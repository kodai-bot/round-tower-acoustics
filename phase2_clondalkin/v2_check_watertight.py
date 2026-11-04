"""
Closed/watertight-manifold check for v2 BEM meshes.

Why this exists: investigation before this rebuild found that v1's
build_tower_mesh.py's classify_skin() only recognises "cyl" and "cone"
face types -- a flat disk cap (e.g. the solid's bottom face at z=0) fails
both tests and is silently discarded (`if kind is None: continue`), so
v1's actual exported .msh is very likely an OPEN tube missing a floor,
not a genuinely closed shell. v1's own validation was done on a sphere,
never on the real tower topology's closedness, so this was never caught.
This is NOT being fixed in v1 (out of scope, v1 stands as-is) -- but v2's
validation gate explicitly requires closed rigid cavities (the analytic
mode targets assume a truly closed-closed cylinder), so v2 must not repeat
the same silent-discard mistake. This script is the explicit check that
catches it.

Method: a closed 2-manifold triangulated surface has every edge shared by
EXACTLY two triangles. Build every triangle's 3 edges as sorted
vertex-index tuples across ALL exported physical groups combined (not
just one domain in isolation -- e.g. a window patch and its neighbouring
wall patch share edges with each other, and both must be counted together
for the whole surface to close up). Any edge appearing once (boundary
edge, i.e. a hole) or more than twice (non-manifold, e.g. duplicated
geometry) is reported as a failure.

Only meshio + numpy, both already in the project .venv -- no new
dependency.
"""
import sys
import numpy as np
import meshio


def is_closed_manifold(mesh_path, verbose=True):
    """Returns (ok: bool, report: dict) for the triangle mesh at mesh_path.

    report contains: n_triangles, n_edges, n_boundary_edges,
    n_nonmanifold_edges, boundary_edge_sample (up to 10 offending edges'
    vertex coordinates, for debugging).
    """
    mesh = meshio.read(mesh_path)
    points = mesh.points

    triangles = []
    for cell_block in mesh.cells:
        if cell_block.type == "triangle":
            triangles.append(cell_block.data)
    if not triangles:
        raise ValueError(f"{mesh_path}: no triangle cells found")
    tris = np.vstack(triangles)

    edge_count = {}
    for tri in tris:
        a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
        for e in ((a, b), (b, c), (c, a)):
            key = (min(e), max(e))
            edge_count[key] = edge_count.get(key, 0) + 1

    boundary_edges = [e for e, n in edge_count.items() if n == 1]
    nonmanifold_edges = [e for e, n in edge_count.items() if n > 2]

    ok = (len(boundary_edges) == 0) and (len(nonmanifold_edges) == 0)

    report = {
        "n_triangles": len(tris),
        "n_edges": len(edge_count),
        "n_boundary_edges": len(boundary_edges),
        "n_nonmanifold_edges": len(nonmanifold_edges),
        "boundary_edge_sample": [
            (points[i].tolist(), points[j].tolist()) for i, j in boundary_edges[:10]
        ],
        "nonmanifold_edge_sample": [
            (points[i].tolist(), points[j].tolist()) for i, j in nonmanifold_edges[:10]
        ],
    }

    if verbose:
        print(f"watertightness check: {mesh_path}")
        print(f"  triangles: {report['n_triangles']}, unique edges: {report['n_edges']}")
        print(f"  boundary (singly-occurring) edges: {report['n_boundary_edges']}")
        print(f"  non-manifold (>2-occurring) edges: {report['n_nonmanifold_edges']}")
        if ok:
            print("  PASS -- closed 2-manifold surface")
        else:
            print("  FAIL -- surface is NOT closed/watertight")
            for p0, p1 in report["boundary_edge_sample"]:
                print(f"    boundary edge: {p0} -- {p1}")
            for p0, p1 in report["nonmanifold_edge_sample"]:
                print(f"    non-manifold edge: {p0} -- {p1}")

    return ok, report


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: v2_check_watertight.py <path-to-.msh>", file=sys.stderr)
        sys.exit(2)
    ok, _ = is_closed_manifold(sys.argv[1])
    sys.exit(0 if ok else 1)
