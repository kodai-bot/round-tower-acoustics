"""
Parametric mesh generator for the v2 (Clondalkin-specific) Round Tower
rebuild -- see v2_params.py for the geometry source and
rebuild_spec.md for the full rationale.

Unlike v1's build_tower_mesh.py (a single zero-thickness shell, windows as
flush Dirichlet patches), v2's wall has REAL thickness (0.745m nominal) and
windows are real rectangular ducts through it. Because of this, interior
air and exterior air are now geometrically distinct volumes joined only
through a duct -- one shared mesh can no longer represent both problems.
This script builds TWO DIFFERENT kinds of mesh depending on --target:

  --target interior : a closed "canister" -- just the r_inner cylinder
      with a real floor AND ceiling (Variant A: rigid flat ceiling at the
      corbel line, cone excluded from the acoustic domain). Used by
      solve_interior_field_v2.py. Independent of wall thickness.
  --target exterior : the real stone exterior -- r_outer cylinder + cone
      roof + floor cap, with the 4 window ducts CUT THROUGH as real holes
      (occ.cut, not occ.fragment) down to r_inner. Used by
      solve_tower_radiation_v2.py. The exterior shape doesn't depend on
      the variant (a hypothetical internal ceiling doesn't change what the
      outdoor air sees), only on wall thickness.

Both meshes use the SAME physical-group convention as v1: domain 1 =
"windows" (the piston/driving boundary), domain 2 = "walls" (rigid/
Neumann). For --target exterior, "windows" is the flat duct-BACK face at
r=r_inner (not the outer mouth) -- see the plan / script docstring in
solve_tower_radiation_v2.py for why the piston is driven at the inner
face: it forces the exterior BEM to actually resolve propagation through
the duct depth, which is the entire point of modelling a real duct
instead of v1's flush patch.

Every generated mesh is checked for closedness (v2_check_watertight.py)
before being trusted -- see that script's docstring for why this check
was added (a real gap found in v1's own mesh, discovered while planning
this rebuild).
"""
import argparse
import math
import sys

import gmsh

import v2_params as P
from v2_check_watertight import is_closed_manifold

# ---------------------------------------------------------------------
# Tolerances for face classification (sampling-grid based, extending v1's
# classify_skin() approach -- see build_tower_mesh.py -- to also catch
# flat disk caps, which v1's version silently discarded).
# ---------------------------------------------------------------------
Z_TOL = 0.05
R_TOL = 0.06
DUCT_LOCAL_TOL = 0.06
FOOTPRINT_ANGLE_SLACK = 1.4   # matches v1's arc_half_width*1.4 margin

OUTER_CUT_MARGIN = 0.3   # m, how far past r_outer the duct-cutting box
                          # extends, purely so the boolean cleanly breaches
                          # the true outer surface (standard OCC practice,
                          # same idea as v1's radial_reach=0.6). The INNER
                          # cut boundary is placed at EXACTLY r_inner (no
                          # margin) -- there's no pre-existing face there
                          # to collide with (the "core" solid is a single
                          # untouched cylinder before cutting), so the
                          # duct-back face lands exactly at r_inner,
                          # matching where the interior solve's driving
                          # pressure is computed.


def window_angles():
    return [i * (2 * math.pi / P.NUM_WINDOWS) for i in range(P.NUM_WINDOWS)]


def sample_face(dim, tag, n=6):
    pmin, pmax = gmsh.model.getParametrizationBounds(dim, tag)
    pts = []
    for i in range(1, n):
        for j in range(1, n):
            u = pmin[0] + (pmax[0] - pmin[0]) * i / n
            v = pmin[1] + (pmax[1] - pmin[1]) * j / n
            pts.append(gmsh.model.getValue(dim, tag, [u, v]))
    return pts


def window_local_coords(x, y, angle):
    """Rotate (x,y) by -angle so the window's outward bearing becomes the
    local +x axis -- i.e. reproduces the box's own pre-rotation frame."""
    ca, sa = math.cos(-angle), math.sin(-angle)
    return x * ca - y * sa, x * sa + y * ca


def nearest_window_hit(x, y, z, angles, sill, head):
    """Returns the window angle this point is plausibly near (angular +
    z footprint, generous slack), or None."""
    r = math.hypot(x, y)
    ang = math.atan2(y, x) % (2 * math.pi)
    if not (sill - Z_TOL <= z <= head + Z_TOL):
        return None
    for wa in angles:
        dang = (ang - wa + math.pi) % (2 * math.pi) - math.pi
        arc_half_width = (P.WINDOW_WIDTH / 2.0) / max(r, 0.3)
        if abs(dang) <= arc_half_width * FOOTPRINT_ANGLE_SLACK:
            return wa
    return None


def classify_interior_face(dim, tag, angles):
    """Interior-A canister: window / floor / ceiling / cyl_inner / None."""
    pts = sample_face(dim, tag)
    if not pts:
        return None
    n = len(pts)
    n_win = n_floor = n_ceil = n_cyl = 0
    for x, y, z in pts:
        r = math.hypot(x, y)
        if abs(z) < Z_TOL:
            n_floor += 1
        if abs(z - P.CHAMBER_HEIGHT) < Z_TOL:
            n_ceil += 1
        if abs(r - P.R_INNER) < R_TOL and -Z_TOL <= z <= P.CHAMBER_HEIGHT + Z_TOL:
            n_cyl += 1
            if nearest_window_hit(x, y, z, angles, P.WINDOW_SILL_Z, P.WINDOW_HEAD_Z) is not None:
                n_win += 1
    if n_win / n > 0.9:
        return "window"
    if n_floor / n > 0.9:
        return "floor"
    if n_ceil / n > 0.9:
        return "ceiling"
    if n_cyl / n > 0.9:
        return "cyl_inner"
    return None


def classify_exterior_face(dim, tag, angles, r_outer, cone_rise):
    """Exterior shell: window (duct-back) / duct_sleeve / floor /
    cyl_outer / cone / None."""
    pts = sample_face(dim, tag)
    if not pts:
        return None
    n = len(pts)
    n_floor = n_outer = n_cone = 0
    n_duct_win = n_duct_sleeve = 0
    for x, y, z in pts:
        r = math.hypot(x, y)
        if abs(z) < Z_TOL:
            n_floor += 1
        if abs(r - r_outer) < R_TOL and -Z_TOL <= z <= P.CHAMBER_HEIGHT + Z_TOL:
            n_outer += 1
        if z > P.CHAMBER_HEIGHT - Z_TOL and cone_rise > 0:
            zc = z - P.CHAMBER_HEIGHT
            r_expected = r_outer * max(0.0, 1.0 - zc / cone_rise)
            if abs(r - r_expected) < R_TOL:
                n_cone += 1
        wa = nearest_window_hit(x, y, z, angles, P.WINDOW_SILL_Z, P.WINDOW_HEAD_Z)
        if wa is not None:
            xl, yl = window_local_coords(x, y, wa)
            if abs(xl - P.R_INNER) < DUCT_LOCAL_TOL:
                n_duct_win += 1
            elif (abs(abs(yl) - P.WINDOW_WIDTH / 2.0) < DUCT_LOCAL_TOL
                  or abs(z - P.WINDOW_SILL_Z) < DUCT_LOCAL_TOL
                  or abs(z - P.WINDOW_HEAD_Z) < DUCT_LOCAL_TOL):
                n_duct_sleeve += 1
    if n_duct_win / n > 0.9:
        return "window"
    if n_duct_sleeve / n > 0.9:
        return "duct_sleeve"
    if n_floor / n > 0.9:
        return "floor"
    if n_outer / n > 0.9:
        return "cyl_outer"
    if n_cone / n > 0.9:
        return "cone"
    return None


def build_interior_mesh(mesh_size, out_path):
    gmsh.initialize()
    gmsh.model.add("v2_tower_interior_A")
    occ = gmsh.model.occ

    plug = occ.addCylinder(0, 0, 0, 0, 0, P.CHAMBER_HEIGHT, P.R_INNER)
    occ.synchronize()

    half_w = P.WINDOW_WIDTH / 2.0
    half_h = P.WINDOW_HEIGHT / 2.0
    radial_reach = 0.6   # box straddles r_inner on both sides, v1-style
    angles = window_angles()
    window_tools = []
    for wa in angles:
        box = occ.addBox(
            P.R_INNER - radial_reach, -half_w, P.WINDOW_SILL_Z,
            2 * radial_reach, P.WINDOW_WIDTH, P.WINDOW_HEIGHT,
        )
        occ.rotate([(3, box)], 0, 0, 0, 0, 0, 1, wa)
        window_tools.append((3, box))
    occ.synchronize()

    pieces, _ = occ.fragment([(3, plug)], window_tools)
    occ.synchronize()

    window_surfs, floor_surfs, ceiling_surfs, wall_surfs = [], [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        kind = classify_interior_face(dim, tag, angles)
        if kind == "window":
            window_surfs.append(tag)
        elif kind == "floor":
            floor_surfs.append(tag)
        elif kind == "ceiling":
            ceiling_surfs.append(tag)
        elif kind == "cyl_inner":
            wall_surfs.append(tag)
        # kind is None -> internal reveal/stub face from the fragment, discard

    all_wall_surfs = floor_surfs + ceiling_surfs + wall_surfs

    print(f"[interior] window surfaces: {window_surfs}")
    print(f"[interior] floor surfaces: {floor_surfs}, ceiling surfaces: {ceiling_surfs}, "
          f"cyl_inner surfaces: {wall_surfs}")
    if not window_surfs:
        print("ERROR: no window surfaces classified", file=sys.stderr)
    if not floor_surfs:
        print("ERROR: no floor surface classified -- cavity would be open!", file=sys.stderr)
    if not ceiling_surfs:
        print("ERROR: no ceiling surface classified -- cavity would be open!", file=sys.stderr)

    gmsh.model.addPhysicalGroup(2, window_surfs, tag=P.WINDOWS_DOMAIN, name="windows")
    gmsh.model.addPhysicalGroup(2, all_wall_surfs, tag=P.WALLS_DOMAIN, name="walls")

    window_area = sum(occ.getMass(2, t) for t in window_surfs)
    wall_area = sum(occ.getMass(2, t) for t in all_wall_surfs)
    expected_window_area = P.NUM_WINDOWS * P.WINDOW_WIDTH * P.WINDOW_HEIGHT
    print(f"[interior] total window area: {window_area:.3f} m^2 (expect ~{expected_window_area:.3f})")
    print(f"[interior] total wall+floor+ceiling area: {wall_area:.3f} m^2")

    _mesh_and_write(mesh_size, out_path)


def build_exterior_mesh(mesh_size, wall_thickness, out_path):
    r_outer = P.R_INNER + wall_thickness
    gmsh.initialize()
    gmsh.model.add("v2_tower_exterior")
    occ = gmsh.model.occ

    solid = occ.addCylinder(0, 0, 0, 0, 0, P.CHAMBER_HEIGHT, r_outer)
    cone = occ.addCone(0, 0, P.CHAMBER_HEIGHT, 0, 0, P.CONE_RISE, r_outer, 0)
    fused, _ = occ.fuse([(3, solid)], [(3, cone)])
    occ.synchronize()

    half_w = P.WINDOW_WIDTH / 2.0
    angles = window_angles()
    duct_tools = []
    radial_span = (r_outer + OUTER_CUT_MARGIN) - P.R_INNER
    for wa in angles:
        box = occ.addBox(
            P.R_INNER, -half_w, P.WINDOW_SILL_Z,
            radial_span, P.WINDOW_WIDTH, P.WINDOW_HEIGHT,
        )
        occ.rotate([(3, box)], 0, 0, 0, 0, 0, 1, wa)
        duct_tools.append((3, box))
    occ.synchronize()

    cut_result, _ = occ.cut(fused, duct_tools)
    occ.synchronize()

    window_surfs, sleeve_surfs, floor_surfs, outer_surfs, cone_surfs = [], [], [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        kind = classify_exterior_face(dim, tag, angles, r_outer, P.CONE_RISE)
        if kind == "window":
            window_surfs.append(tag)
        elif kind == "duct_sleeve":
            sleeve_surfs.append(tag)
        elif kind == "floor":
            floor_surfs.append(tag)
        elif kind == "cyl_outer":
            outer_surfs.append(tag)
        elif kind == "cone":
            cone_surfs.append(tag)
        # None -> discard

    all_wall_surfs = sleeve_surfs + floor_surfs + outer_surfs + cone_surfs

    print(f"[exterior] window (duct-back) surfaces: {window_surfs}")
    print(f"[exterior] duct_sleeve: {sleeve_surfs}, floor: {floor_surfs}, "
          f"cyl_outer: {outer_surfs}, cone: {cone_surfs}")
    if not window_surfs:
        print("ERROR: no window (duct-back) surfaces classified", file=sys.stderr)
    if len(sleeve_surfs) < P.NUM_WINDOWS:
        print(f"WARNING: expected >= {P.NUM_WINDOWS} duct_sleeve faces (4 per duct "
              f"typically: 2 sides + top + bottom), got {len(sleeve_surfs)}", file=sys.stderr)
    if not floor_surfs:
        print("ERROR: no floor surface classified -- shell would be open!", file=sys.stderr)

    gmsh.model.addPhysicalGroup(2, window_surfs, tag=P.WINDOWS_DOMAIN, name="windows")
    gmsh.model.addPhysicalGroup(2, all_wall_surfs, tag=P.WALLS_DOMAIN, name="walls")

    window_area = sum(occ.getMass(2, t) for t in window_surfs)
    wall_area = sum(occ.getMass(2, t) for t in all_wall_surfs)
    expected_window_area = P.NUM_WINDOWS * P.WINDOW_WIDTH * P.WINDOW_HEIGHT
    print(f"[exterior] total window (duct-back) area: {window_area:.3f} m^2 "
          f"(expect ~{expected_window_area:.3f})")
    print(f"[exterior] total wall+floor+cone+sleeve area: {wall_area:.3f} m^2")

    _mesh_and_write(mesh_size, out_path)


def _mesh_and_write(mesh_size, out_path):
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size / 3)
    gmsh.model.mesh.generate(2)

    n_elems_by_group = {}
    for dim, tag in gmsh.model.getPhysicalGroups(2):
        ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        count = 0
        for e in ents:
            elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim, e)
            for et in elem_tags:
                count += len(et)
        n_elems_by_group[tag] = count
    print(f"elements per physical group: {n_elems_by_group}")

    gmsh.write(out_path)
    print(f"wrote {out_path}")
    gmsh.finalize()

    ok, report = is_closed_manifold(out_path)
    if not ok:
        print(f"FATAL: {out_path} is not watertight -- refusing to consider this mesh usable. "
              f"{report['n_boundary_edges']} boundary edges, "
              f"{report['n_nonmanifold_edges']} non-manifold edges.", file=sys.stderr)
        sys.exit(1)
    print(f"watertightness check passed for {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", required=True, choices=["interior", "exterior"])
    ap.add_argument("--variant", default=P.VARIANT, choices=["A"],
                     help="Only Variant A (ceiled) is implemented this round; "
                          "Variant B is deferred to a follow-up round.")
    ap.add_argument("--wall-thickness", type=float, default=P.WALL_THICKNESS,
                     help="Exterior mesh only. Sensitivity cases (0.64/0.813m) "
                          "are deferred this round but this flag is wired up "
                          "for that follow-up work.")
    ap.add_argument("--mesh-size", type=float, default=P.MESH_SIZE)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.target == "interior":
        out_path = args.out or P.MESH_INTERIOR_A_PATH
        build_interior_mesh(args.mesh_size, out_path)
    else:
        out_path = args.out or P.MESH_EXTERIOR_PATH
        build_exterior_mesh(args.mesh_size, args.wall_thickness, out_path)


if __name__ == "__main__":
    main()
