"""
Parametric mesh generator for the Round Tower upper chamber, for exterior
acoustic radiation (BEM) analysis.

Model philosophy (matches Computational_Towers.odt methodology):
  - We only need the CLOSED exterior skin of the chamber as a single thin
    surface (cylinder wall + roof cap, now a cone -- see below). Wall
    thickness and the interior air volume are irrelevant for the exterior-
    radiation BEM problem, so we do NOT model a hollow wall or a floor --
    just the outer silhouette.
  - Windows are NOT real holes. Following the paper's own formulation, they
    are Dirichlet-condition patches (prescribed pressure from the interior
    source) embedded flush in an otherwise closed, sound-hard (Neumann)
    shell. So the geometry step just needs to PARTITION the skin into two
    physical surface groups -- "windows" and "walls" -- with no material
    removed and no true apertures cut. This avoids the open-screen BEM
    problem entirely and matches tower_with_physical_groups.msh's original
    (if buggy) intent.

Dimensions recovered from the embedded generator script inside
tower_room_2.blend (see project chat log). Values are deliberately rough --
survey dimensions for these towers are not reliable, so this is a general
model, not a precise reconstruction. Adjust the constants below freely.
"""

import os
HERE = os.path.dirname(os.path.abspath(__file__))

import gmsh
import math
import sys

# ---------------------------------------------------------------------
# Parameters (all metres / radians unless noted)
# ---------------------------------------------------------------------
TOWER_RADIUS   = 2.2    # representative single-surface radius (recovered
                         # model had inner 1.9 / outer 2.5; BEM shell uses
                         # one thin surface, so we split the difference).
                         # No new site-specific radius data yet -- unchanged.

# --- Roof geometry: REVISED from a shallow-dome guess to a CONE, per an
# actual (imperfect) 3D scan of Clondalkin round tower's intact interior
# roof: "a regular cone of smooth corbels," springing (cone base) sitting
# just above the window tops (essentially no blank wall between window top
# and roof start), rising ~8ft above the window tops to a point apex.
# This replaces the earlier dome assumption entirely for this run.
WALL_HEIGHT    = 1.2192  # = WINDOW_HEIGHT below -- springing sits right at
                          # window-top height ("cone begins to slope just
                          # above the windows", read here as ~zero margin)
CONE_RISE      = 2.5908  # 8ft6in, apex above the springing/window-top
                          # (the scan gave two mutually-consistent figures
                          # for this: "8ft6in down to the top of the window
                          # opening" measured from the apex, and "over 8
                          # feet above the top of the windows" restated
                          # independently -- used the more specific one)

NUM_WINDOWS    = 4
WINDOW_WIDTH   = 0.6096  # 2ft, Clondalkin scan ("about 2 foot across") --
                         # supersedes the earlier 0.46m general-survey
                         # average for this run; site-specific, but from
                         # an explicitly imperfect scan
WINDOW_HEIGHT  = 1.2192  # 4ft, Clondalkin scan ("about 4 foot high",
                         # stated directly). NOTE: the scan's two
                         # apex-referenced measurements (8ft6in to window
                         # top, 12ft to window bottom) subtract to 3ft6in
                         # (1.067m), ~0.15m shorter than the direct 4ft
                         # estimate -- a real discrepancy, consistent with
                         # "imperfect 3D scan," not silently resolved. Used
                         # the direct estimate; if the subtracted figure is
                         # more trustworthy this should be revisited.
WINDOW_BASE_Z  = 0.0    # base at floor level (unchanged assumption --
                         # this scan didn't re-confirm floor-level bases,
                         # only roof/window-top geometry)

MESH_SIZE      = 0.15   # ~15cm elements -> ~5.7 elements/wavelength at
                         # 400Hz. Originally tried 0.1m (~8.6 elem/wl at
                         # 400Hz, comfortably within the paper's 6-10 rule)
                         # but that pushed element count to ~13,500 and
                         # dense BEM assembly cost (no FMM/OpenCL available
                         # in this environment) blew up to ~300s PER
                         # FREQUENCY SOLVE -- confirmed via direct timing,
                         # not GMRES-tolerance-dependent, purely assembly
                         # cost scaling with element count. 0.15m is a
                         # deliberate compromise: ~50s/solve, same margin
                         # (just under 6 elem/wavelength) at its own 400Hz
                         # ceiling that the original 0.2m mesh had at 300Hz.

OUT_DIR = f"{HERE}"

# ---------------------------------------------------------------------
gmsh.initialize()
gmsh.model.add("round_tower_upper_chamber")
occ = gmsh.model.occ

# --- Build the exterior silhouette as a solid "plug" -------------------
# We don't care about the interior of this solid at all -- it exists only
# so OpenCASCADE's boolean kernel can robustly hand us a clean, correctly
# split outer skin. No hollowing, no wall thickness, no floor.

cyl = occ.addCylinder(0, 0, 0, 0, 0, WALL_HEIGHT, TOWER_RADIUS)

# Cone roof: a straight-sided cone from radius TOWER_RADIUS at the
# springing (z=WALL_HEIGHT) tapering linearly to a point apex at
# z=WALL_HEIGHT+CONE_RISE. addCone(x,y,z, dx,dy,dz, r1,r2) -- r2=0 gives a
# true point apex, matching "smooth corbels up to the apex."
cone = occ.addCone(0, 0, WALL_HEIGHT, 0, 0, CONE_RISE, TOWER_RADIUS, 0)

tower, _ = occ.fuse([(3, cyl)], [(3, cone)])
occ.synchronize()

# --- Window cutters: boxes that straddle the skin at each cardinal point
# so a fragment splits off exactly the right patch of curved skin.
half_w = WINDOW_WIDTH / 2.0
half_h = WINDOW_HEIGHT / 2.0
radial_reach = 0.6  # box just needs to straddle r=R on both sides so the
                     # fragment cleanly slices through the skin
window_tools = []
for i in range(NUM_WINDOWS):
    angle = i * (2 * math.pi / NUM_WINDOWS)
    box = occ.addBox(
        TOWER_RADIUS - radial_reach, -half_w, WINDOW_BASE_Z,
        2 * radial_reach, WINDOW_WIDTH, WINDOW_HEIGHT,
    )
    occ.rotate([(3, box)], 0, 0, 0, 0, 0, 1, angle)
    window_tools.append((3, box))

occ.synchronize()

all_pieces, _ = occ.fragment(tower, window_tools)
occ.synchronize()

# ---------------------------------------------------------------------
# Classify resulting surfaces. Fragmenting a solid with boxes that punch
# through it creates a lot of geometry we don't want: internal "reveal"
# faces (the sides of each box's intersection with the tower body) whose
# radius sweeps from the interior out to the skin, plus stub faces from
# the parts of the boxes sticking out past the skin into open air. Only
# a face that is a genuine (thin) piece of the original cylinder/dome
# skin -- i.e. one whose radius stays close to TOWER_RADIUS (or lies on
# the dome ellipsoid) EVERYWHERE across it, not just at its centroid --
# is real skin. We test this by sampling each face's parametrisation on
# a grid rather than trusting its centroid, since the centroid of a
# wraparound face (like the main wall band, which still touches all 360
# degrees around the 4 window notches) sits near the central axis and is
# useless as a classifier on its own.
# ---------------------------------------------------------------------

def in_window_footprint(x, y, z):
    r = math.hypot(x, y)
    ang = math.atan2(y, x) % (2 * math.pi)
    for i in range(NUM_WINDOWS):
        wangle = i * (2 * math.pi / NUM_WINDOWS)
        dang = (ang - wangle + math.pi) % (2 * math.pi) - math.pi
        arc_half_width = half_w / TOWER_RADIUS
        if abs(dang) <= arc_half_width * 1.4 and WINDOW_BASE_Z - 0.05 <= z <= WINDOW_BASE_Z + WINDOW_HEIGHT + 0.05:
            if abs(r - TOWER_RADIUS) < 0.4:
                return True
    return False


def sample_face(dim, tag, n=6):
    pmin, pmax = gmsh.model.getParametrizationBounds(dim, tag)
    pts = []
    for i in range(1, n):
        for j in range(1, n):
            u = pmin[0] + (pmax[0] - pmin[0]) * i / n
            v = pmin[1] + (pmax[1] - pmin[1]) * j / n
            pts.append(gmsh.model.getValue(dim, tag, [u, v]))
    return pts


def classify_skin(dim, tag):
    """Returns 'cyl', 'cone', or None (not a real skin patch)."""
    pts = sample_face(dim, tag)
    if not pts:
        return None
    n_cyl = n_cone = 0
    for x, y, z in pts:
        r = math.hypot(x, y)
        if abs(r - TOWER_RADIUS) < 0.08 and -0.05 <= z <= WALL_HEIGHT + 0.05:
            n_cyl += 1
        if z > WALL_HEIGHT - 0.08 and CONE_RISE > 0:
            zc = z - WALL_HEIGHT
            r_expected = TOWER_RADIUS * max(0.0, 1.0 - zc / CONE_RISE)
            if abs(r - r_expected) < 0.08:
                n_cone += 1
    if n_cyl / len(pts) > 0.9:
        return "cyl"
    if n_cone / len(pts) > 0.9:
        return "cone"
    return None


window_surfs = []
wall_surfs = []

for dim, tag in gmsh.model.getEntities(2):
    kind = classify_skin(dim, tag)
    if kind is None:
        continue  # internal reveal face or exterior stub -- discard
    com = gmsh.model.occ.getCenterOfMass(dim, tag)
    if kind == "cyl" and in_window_footprint(*com):
        window_surfs.append(tag)
    else:
        wall_surfs.append(tag)

print(f"window surfaces: {window_surfs}")
print(f"wall surfaces:   {wall_surfs}")

if len(window_surfs) == 0:
    print("ERROR: no window surfaces classified -- geometry/logic problem", file=sys.stderr)

windows_pg = gmsh.model.addPhysicalGroup(2, window_surfs, name="windows")
walls_pg = gmsh.model.addPhysicalGroup(2, wall_surfs, name="walls")

window_area = sum(gmsh.model.occ.getMass(2, t) for t in window_surfs)
wall_area = sum(gmsh.model.occ.getMass(2, t) for t in wall_surfs)
print(f"total window area: {window_area:.3f} m^2 (expect ~= "
      f"{NUM_WINDOWS * WINDOW_WIDTH * WINDOW_HEIGHT:.3f} m^2)")
print(f"total wall+dome area: {wall_area:.3f} m^2")

# ---------------------------------------------------------------------
# Mesh only these skin surfaces (2D surface mesh -- BEM needs no volume
# mesh at all). Mesh.SaveAll=0 means only elements on surfaces that
# belong to a physical group get written out, so the leftover internal
# reveal/stub faces from the fragment never make it into the output file
# even though they still exist (and get meshed) in memory.
# ---------------------------------------------------------------------
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.option.setNumber("Mesh.MeshSizeMax", MESH_SIZE)
gmsh.option.setNumber("Mesh.MeshSizeMin", MESH_SIZE / 3)
gmsh.model.mesh.generate(2)

n_win_elems = sum(len(gmsh.model.mesh.getElements(2, t)[1][0])
                   if gmsh.model.mesh.getElements(2, t)[1] else 0
                   for t in window_surfs)
n_wall_elems = sum(len(gmsh.model.mesh.getElements(2, t)[1][0])
                    if gmsh.model.mesh.getElements(2, t)[1] else 0
                    for t in wall_surfs)
print(f"window elements: {n_win_elems}, wall elements: {n_wall_elems}")

out_path = f"{OUT_DIR}/tower_upper_chamber.msh"
gmsh.write(out_path)
print(f"wrote {out_path}")

gmsh.finalize()
