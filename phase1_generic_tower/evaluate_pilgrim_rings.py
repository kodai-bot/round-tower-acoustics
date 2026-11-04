"""
Evaluate exterior SPL at points a pilgrim would actually walk through,
circumambulating the 2nd and 3rd concentric rings of the monastic
enclosure -- rather than the generic tower-centred distance/azimuth grid
used in solve_tower_radiation.py.

Site geometry (from the project description):
  - Common central origin = centre of the enclosure (tabernacle/church).
  - 1st ring (clergy/inner sanctum) radius 20m -- not evaluated here,
    pilgrims aren't in this zone.
  - 2nd ring (pilgrims, prayer mats, circumambulation) radius 45m.
  - 3rd ring (outer boundary) radius 91.5m.
  - Tower sits at the INNER EDGE of the 2nd ring, WEST of the church
    doorway, 17m from the common origin. Tower door (and, by the same
    facing logic, window 0 in our BEM model) faces EAST toward the
    centre of the complex.

Because the tower is offset from the rings' common centre (not
concentric with them), distance from the tower to a pilgrim varies
continuously as they walk around a ring -- from ~28m (2nd ring, tower
side) to ~62m (2nd ring, far side), and ~74.5m to ~108.5m for the 3rd
ring. This script walks each ring point-by-point in the SITE frame, and
transforms into the tower's own local frame for BEM evaluation.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential

MESH_PATH = f"{HERE}/tower_upper_chamber.msh"
DRIVING_CSV = f"{HERE}/window_driving_conditions.csv"
OUT_DIR = f"{HERE}"

c = 343.0
rho0 = 1.225
WINDOWS_DOMAIN = 1
NUM_WINDOWS = 4
WINDOW_REF_ANGLES = np.array([i * 360 / NUM_WINDOWS for i in range(NUM_WINDOWS)])

TOWER_OFFSET = 17.0     # metres, tower centre from the common origin, due west
RING_RADII = {"2nd_ring": 45.0, "3rd_ring": 91.5}
N_POINTS_AROUND_RING = 36  # every 10 degrees

TOWER_ELEVATION = 30.0
LISTENER_HEIGHT = 1.6

FREQUENCIES = [50, 75, 100, 125, 143, 150, 160, 165, 170, 175, 180, 185, 190,
               196, 200, 207, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
               310, 320, 330, 340, 349, 360, 375, 390]
SOURCE_CASES = ["standing"]  # "seated" control not available at this mesh/
                              # frequency set -- see solve_interior_field.py

# ---------------------------------------------------------------------
grid = bempp.import_grid(MESH_PATH)
space = bempp.function_space(grid, "DP", 0)
identity = bempp.operators.boundary.sparse.identity(space, space, space)
elem_volumes = grid.volumes
window_elems = grid.domain_indices == WINDOWS_DOMAIN
window_area_each = float(elem_volumes[window_elems].sum()) / NUM_WINDOWS

driving = {}
with open(DRIVING_CSV) as fh:
    for row in csv.DictReader(fh):
        key = (row["source_case"], int(float(row["freq_hz"])))
        driving.setdefault(key, [None] * NUM_WINDOWS)
        driving[key][int(row["window_index"])] = complex(float(row["pressure_re"]), float(row["pressure_im"]))


def solve_at_frequency(f, source_case):
    omega = 2 * np.pi * f
    k = omega / c
    dlp = bempp.operators.boundary.helmholtz.double_layer(space, space, space, k)
    slp = bempp.operators.boundary.helmholtz.single_layer(space, space, space, k)
    lhs = 0.5 * identity - dlp

    p_windows = np.array(driving[(source_case, f)], dtype=np.complex128)
    ref_angles = WINDOW_REF_ANGLES.copy()

    @bempp.complex_callable
    def neumann_data(x, n, domain_index, result):
        if domain_index != WINDOWS_DOMAIN:
            result[0] = 0.0
            return
        az = np.degrees(np.arctan2(x[1], x[0])) % 360.0
        best_w, best_d = 0, 1e9
        for wi in range(ref_angles.shape[0]):
            d = np.abs((az - ref_angles[wi] + 180.0) % 360.0 - 180.0)
            if d < best_d:
                best_d, best_w = d, wi
        result[0] = 1j * k * p_windows[best_w]

    g_fun = bempp.GridFunction(space, fun=neumann_data)
    rhs = -1.0 * (slp * g_fun)
    p_sol, info = gmres(lhs, rhs, tol=1e-5)
    if info != 0:
        print(f"  [WARNING] f={f}Hz source={source_case}: gmres info={info}")

    Q = sum(abs(p) / (rho0 * c) * window_area_each for p in p_windows)
    return k, p_sol, g_fun, Q


def spl_db(p_complex):
    p_rms = np.abs(p_complex) / np.sqrt(2)
    return 20 * np.log10(np.maximum(p_rms, 1e-12) / 20e-6)


# ---------------------------------------------------------------------
for source_case in SOURCE_CASES:
    rows = []
    for f in FREQUENCIES:
        if (source_case, f) not in driving:
            continue
        print(f"solving f={f}Hz source={source_case} ...")
        k, p_sol, g_fun, Q = solve_at_frequency(f, source_case)

        for ring_name, R in RING_RADII.items():
            thetas = np.linspace(0, 360, N_POINTS_AROUND_RING, endpoint=False)
            # site-frame ring point -> tower-local frame (tower sits at
            # (-TOWER_OFFSET, 0) in site coords, window 0 / door faces
            # +x site-east, i.e. tower-local axes align with site axes)
            site_x = R * np.cos(np.radians(thetas))
            site_y = R * np.sin(np.radians(thetas))
            local_x = site_x + TOWER_OFFSET
            local_y = site_y
            local_z = np.full_like(local_x, -TOWER_ELEVATION + LISTENER_HEIGHT)
            dist_from_tower = np.sqrt(local_x ** 2 + local_y ** 2)

            pts = np.vstack([local_x, local_y, local_z])
            dlp_pot = helmholtz_potential.double_layer(space, pts, k)
            slp_pot = helmholtz_potential.single_layer(space, pts, k)
            p_tower = (dlp_pot.evaluate(p_sol) - slp_pot.evaluate(g_fun)).flatten()

            r_slant = np.sqrt(local_x ** 2 + local_y ** 2 + local_z ** 2)
            omega = k * c
            p_free = 1j * omega * rho0 * Q / (4 * np.pi * r_slant) * np.exp(1j * k * r_slant)

            for th, d, pt, pf in zip(thetas, dist_from_tower, p_tower, p_free):
                rows.append({
                    "freq_hz": f, "ring": ring_name, "theta_site_deg": th,
                    "dist_from_tower_m": d,
                    "spl_tower_db": spl_db(pt), "spl_free_db": spl_db(pf),
                    "gain_db": spl_db(pt) - spl_db(pf),
                })

    out_csv = f"{OUT_DIR}/pilgrim_ring_walk_{source_case}.csv"
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_csv} ({len(rows)} rows)")
