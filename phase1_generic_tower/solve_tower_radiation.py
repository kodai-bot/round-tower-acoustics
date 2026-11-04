"""
Exterior acoustic radiation solver for the Round Tower upper chamber.

Physics / modelling choices (read this before trusting the numbers):

  - Boundary integral equation solved (direct method, non-combined):
        (0.5*Id - K) p = -V g
    where p is the (unknown) surface pressure, g = i*omega*rho0*v_n is
    the prescribed Neumann data, K = double layer, V = single layer.
    Validated against the exact pulsating-sphere solution to <0.6% error
    (see validate_bem_formulation.py) before being trusted here.
    CAVEAT: this direct (non-CFIE) formulation can fail at isolated
    "irregular frequencies" coinciding with interior Dirichlet eigenmodes
    of the volume the closed shell encloses. We flag (not silently trust)
    any frequency where GMRES doesn't converge cleanly.

  - Field at exterior points x is then p(x) = D[p](x) - S[g](x).

  - Window driving data now comes from solve_interior_field.py -- a
    validated interior rigid-cavity BEM solve driven by an omnidirectional
    point-source monopole (standing in for a chanting voice), read off
    individually at each of the 4 windows. This replaces the old flat
    guessed-SPL approach entirely: each window gets its own computed
    complex pressure (magnitude AND phase), varying by frequency and by
    which of the two source cases (standing-central / seated-control) is
    being tested. Converted to a velocity boundary condition via the same
    plane-wave impedance relation used throughout this project
    (v = p/(rho0*c)), giving g = i*k*p_window at each window.

  - This is a WEAKLY COUPLED, two-step approximation: the interior field
    was computed treating the windows as closed/rigid (ignoring the ~5%
    area they actually occupy), then that pressure drives the exterior
    problem as if the window were a vibrating piston. It does not solve
    a single self-consistent interior-exterior transmission problem.
    Reasonable given the window area fraction, but not exact.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential

# ---------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------
MESH_PATH = f"{HERE}/tower_upper_chamber.msh"
DRIVING_CSV = f"{HERE}/window_driving_conditions.csv"
OUT_DIR = f"{HERE}"

c = 343.0        # speed of sound, m/s
rho0 = 1.225     # air density, kg/m^3

WINDOWS_DOMAIN = 1
WALLS_DOMAIN = 2
NUM_WINDOWS = 4
WINDOW_REF_ANGLES = np.array([i * 360 / NUM_WINDOWS for i in range(NUM_WINDOWS)])

TOWER_ELEVATION = 30.0  # metres, upper chamber floor above ground (your
                         # 90-110ft estimate -> ~27-34m; 30m is the middle)
LISTENER_HEIGHT = 1.6    # metres, head height of a listener on the ground

FREQUENCIES = [50, 75, 100, 125, 143, 150, 160, 165, 170, 175, 180, 185, 190,
               196, 200, 207, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
               310, 320, 330, 340, 349, 360, 375, 390]  # Hz -- extended to
# cover French's 5th/6th scale degrees; mesh re-cut to 0.15m elements
# NOTE: the mesh was generated at 0.2m elements (~6-10 elem/wavelength up
# to 300Hz per Computational_Towers.odt's own target). Frequencies above
# 300Hz are NOT resolved by this mesh -- the 2.5-3.5kHz singer-formant
# band from the paper needs a much finer mesh and is deliberately
# excluded here.

GROUND_DISTANCES = [50, 100, 200, 500, 1000, 2000]  # metres from tower base
AZIMUTHS_DEG = list(range(0, 360, 10))  # 36 points, 10 deg spacing --
# widened from an earlier 8-point/45deg grid: the tower has exact 4-fold
# (90deg) window symmetry, and 45deg sampling is a submultiple of that
# symmetry, risking aliasing with the real lobe structure. Confirmed real:
# at 250Hz/100m the true pattern alternates between ~-20.7dB exactly on a
# window axis (0/90/180/270deg) and ~-3.6dB on the diagonals -- a ~17dB
# swing an 8-point grid can (and did) capture at the SPECIFIC points it
# happened to sample, but finer resolution is needed to confirm there's no
# additional structure between those points. 0 deg aligned with window 0's
# +x axis.

ELEVATION_RADIUS = 200.0  # metres, fixed slant range for the
                           # horizontal/downward vs upward directivity check
ELEVATION_ANGLES_DEG = list(range(-90, 91, 15))  # -90 = straight down,
                                                   # 0 = horizontal,
                                                   # +90 = straight up
ELEVATION_AZIMUTHS_DEG = [0, 45, 90, 135, 180, 225, 270, 315]  # kept at the
# original 8-point resolution for the elevation sweep specifically (not
# widened -- Task 3 only asked for finer azimuth resolution in the ground-
# projection sweep; widening this too would roughly 4.5x the elevation
# sweep's point count for no requested benefit)

# Pilgrim ring-walk geometry (site layout, see evaluate_pilgrim_rings.py
# for the original standalone version -- merged in here so the ring-walk
# evaluation reuses the SAME solved p_sol/g_fun per frequency instead of
# re-solving the whole BEM problem a second time. That earlier standalone
# script cost ~25-30 minutes per run purely by duplicating this solve for
# different evaluation points; evaluating extra points from an
# already-solved field is nearly free by comparison.)
TOWER_OFFSET = 17.0
RING_RADII = {"2nd_ring": 45.0, "3rd_ring": 91.5}
N_POINTS_AROUND_RING = 36

SOURCE_CASES = ["standing"]  # "seated" control skipped for this extended-
# range (0.15m mesh) run -- see note in solve_interior_field.py

# ---------------------------------------------------------------------
grid = bempp.import_grid(MESH_PATH)
space = bempp.function_space(grid, "DP", 0)
identity = bempp.operators.boundary.sparse.identity(space, space, space)

window_elems = grid.domain_indices == WINDOWS_DOMAIN
wall_elems = grid.domain_indices == WALLS_DOMAIN
print(f"mesh: {grid.number_of_elements} elements "
      f"({window_elems.sum()} window, {wall_elems.sum()} wall/dome)")
elem_volumes = grid.volumes
window_area_total = float(elem_volumes[window_elems].sum())
window_area_each = window_area_total / NUM_WINDOWS

# load computed per-window driving pressures: (source_case, freq) -> [p0..p3]
driving = {}
with open(DRIVING_CSV) as fh:
    for row in csv.DictReader(fh):
        key = (row["source_case"], int(float(row["freq_hz"])))
        driving.setdefault(key, [None] * NUM_WINDOWS)
        w = int(row["window_index"])
        driving[key][w] = complex(float(row["pressure_re"]), float(row["pressure_im"]))


def nearest_window_index(az_deg):
    diffs = np.abs((az_deg - WINDOW_REF_ANGLES + 180) % 360 - 180)
    return int(np.argmin(diffs))


def solve_at_frequency(f, source_case):
    omega = 2 * np.pi * f
    k = omega / c

    dlp = bempp.operators.boundary.helmholtz.double_layer(space, space, space, k)
    slp = bempp.operators.boundary.helmholtz.single_layer(space, space, space, k)
    lhs = 0.5 * identity - dlp

    # numba (used by bempp's @complex_callable JIT) can't call back into
    # arbitrary Python functions, so the window-bucketing logic has to be
    # inlined directly rather than calling nearest_window_index().
    p_windows = np.array(driving[(source_case, f)], dtype=np.complex128)
    ref_angles = WINDOW_REF_ANGLES.copy()

    @bempp.complex_callable
    def neumann_data(x, n, domain_index, result):
        if domain_index != WINDOWS_DOMAIN:
            result[0] = 0.0
            return
        az = np.degrees(np.arctan2(x[1], x[0])) % 360.0
        best_w = 0
        best_d = 1e9
        for wi in range(ref_angles.shape[0]):
            d = np.abs((az - ref_angles[wi] + 180.0) % 360.0 - 180.0)
            if d < best_d:
                best_d = d
                best_w = wi
        result[0] = 1j * k * p_windows[best_w]

    g_fun = bempp.GridFunction(space, fun=neumann_data)
    rhs = -1.0 * (slp * g_fun)

    p_sol, info = gmres(lhs, rhs, tol=1e-5)
    if info != 0:
        print(f"  [WARNING] f={f}Hz source={source_case}: gmres did not "
              f"cleanly converge (info={info}) -- possible irregular "
              f"frequency, treat result with suspicion")

    # fair free-field baseline: a monopole with the same total volume
    # velocity actually flowing out of the windows in this solve
    Q = sum(abs(p) / (rho0 * c) * window_area_each for p in p_windows)

    return k, omega, p_sol, g_fun, Q


def evaluate_ground_points(k, p_sol, g_fun, Q):
    pts = []
    meta = []
    for d in GROUND_DISTANCES:
        for az_deg in AZIMUTHS_DEG:
            az = np.radians(az_deg)
            x = d * np.cos(az)
            y = d * np.sin(az)
            z = -TOWER_ELEVATION + LISTENER_HEIGHT
            pts.append((x, y, z))
            meta.append((d, az_deg))
    pts = np.array(pts).T

    dlp_pot = helmholtz_potential.double_layer(space, pts, k)
    slp_pot = helmholtz_potential.single_layer(space, pts, k)
    p_tower = (dlp_pot.evaluate(p_sol) - slp_pot.evaluate(g_fun)).flatten()

    r_slant = np.sqrt(pts[0] ** 2 + pts[1] ** 2 + pts[2] ** 2)
    omega = k * c
    p_free = 1j * omega * rho0 * Q / (4 * np.pi * r_slant) * np.exp(1j * k * r_slant)

    return meta, p_tower, p_free


def evaluate_elevation_sweep(k, p_sol, g_fun, Q):
    pts = []
    meta = []
    for az_deg in ELEVATION_AZIMUTHS_DEG:
        for el_deg in ELEVATION_ANGLES_DEG:
            az = np.radians(az_deg)
            el = np.radians(el_deg)
            horiz = ELEVATION_RADIUS * np.cos(el)
            x = horiz * np.cos(az)
            y = horiz * np.sin(az)
            z = ELEVATION_RADIUS * np.sin(el)  # relative to the tower's own
                                                 # floor level (BEM origin) --
                                                 # tower height (~3.4m) is
                                                 # negligible against 200m
            pts.append((x, y, z))
            meta.append((az_deg, el_deg))
    pts = np.array(pts).T

    dlp_pot = helmholtz_potential.double_layer(space, pts, k)
    slp_pot = helmholtz_potential.single_layer(space, pts, k)
    p_tower = (dlp_pot.evaluate(p_sol) - slp_pot.evaluate(g_fun)).flatten()

    omega = k * c
    p_free = 1j * omega * rho0 * Q / (4 * np.pi * ELEVATION_RADIUS) * np.exp(1j * k * ELEVATION_RADIUS)
    p_free_arr = np.full_like(p_tower, p_free)

    return meta, p_tower, p_free_arr


def evaluate_pilgrim_rings(k, p_sol, g_fun, Q):
    pts = []
    meta = []
    for ring_name, R in RING_RADII.items():
        thetas = np.linspace(0, 360, N_POINTS_AROUND_RING, endpoint=False)
        site_x = R * np.cos(np.radians(thetas))
        site_y = R * np.sin(np.radians(thetas))
        local_x = site_x + TOWER_OFFSET
        local_y = site_y
        local_z = np.full_like(local_x, -TOWER_ELEVATION + LISTENER_HEIGHT)
        dist_from_tower = np.sqrt(local_x ** 2 + local_y ** 2)
        for th, d, lx, ly, lz in zip(thetas, dist_from_tower, local_x, local_y, local_z):
            pts.append((lx, ly, lz))
            meta.append((ring_name, th, d))
    pts = np.array(pts).T

    dlp_pot = helmholtz_potential.double_layer(space, pts, k)
    slp_pot = helmholtz_potential.single_layer(space, pts, k)
    p_tower = (dlp_pot.evaluate(p_sol) - slp_pot.evaluate(g_fun)).flatten()

    r_slant = np.sqrt(pts[0] ** 2 + pts[1] ** 2 + pts[2] ** 2)
    omega = k * c
    p_free = 1j * omega * rho0 * Q / (4 * np.pi * r_slant) * np.exp(1j * k * r_slant)

    return meta, p_tower, p_free


def spl_db(p_complex):
    p_rms = np.abs(p_complex) / np.sqrt(2)
    return 20 * np.log10(np.maximum(p_rms, 1e-12) / 20e-6)


def intensity(p_complex):
    """proportional-to-intensity weight for a crude hemisphere power split"""
    return np.abs(p_complex) ** 2


# ---------------------------------------------------------------------
for source_case in SOURCE_CASES:
    ground_rows = []
    elev_rows = []
    ring_rows = []
    for f in FREQUENCIES:
        if (source_case, f) not in driving:
            print(f"skipping f={f} source={source_case}: no driving data")
            continue
        print(f"solving f={f} Hz source={source_case} ...")
        k, omega, p_sol, g_fun, Q = solve_at_frequency(f, source_case)

        meta, p_tower, p_free = evaluate_ground_points(k, p_sol, g_fun, Q)
        for (d, az), pt, pf in zip(meta, p_tower, p_free):
            ground_rows.append({
                "freq_hz": f, "distance_m": d, "azimuth_deg": az,
                "spl_tower_db": spl_db(pt), "spl_free_db": spl_db(pf),
                "gain_db": spl_db(pt) - spl_db(pf),
            })

        emeta, ep_tower, ep_free = evaluate_elevation_sweep(k, p_sol, g_fun, Q)
        for (az, el), pt, pf in zip(emeta, ep_tower, ep_free):
            elev_rows.append({
                "freq_hz": f, "azimuth_deg": az, "elevation_deg": el,
                "spl_tower_db": spl_db(pt), "spl_free_db": spl_db(pf),
                "gain_db": spl_db(pt) - spl_db(pf),
            })

        rmeta, rp_tower, rp_free = evaluate_pilgrim_rings(k, p_sol, g_fun, Q)
        for (ring, theta, dist), pt, pf in zip(rmeta, rp_tower, rp_free):
            ring_rows.append({
                "freq_hz": f, "ring": ring, "theta_site_deg": theta,
                "dist_from_tower_m": dist,
                "spl_tower_db": spl_db(pt), "spl_free_db": spl_db(pf),
                "gain_db": spl_db(pt) - spl_db(pf),
            })

    out_csv = f"{OUT_DIR}/ground_projection_{source_case}.csv"
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(ground_rows[0].keys()))
        writer.writeheader()
        writer.writerows(ground_rows)
    print(f"wrote {out_csv} ({len(ground_rows)} rows)")

    out_csv2 = f"{OUT_DIR}/elevation_sweep_{source_case}.csv"
    with open(out_csv2, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(elev_rows[0].keys()))
        writer.writeheader()
        writer.writerows(elev_rows)
    print(f"wrote {out_csv2} ({len(elev_rows)} rows)")

    out_csv3 = f"{OUT_DIR}/pilgrim_ring_walk_{source_case}.csv"
    with open(out_csv3, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(ring_rows[0].keys()))
        writer.writeheader()
        writer.writerows(ring_rows)
    print(f"wrote {out_csv3} ({len(ring_rows)} rows)")

    # rough hemisphere power split (equal-angle grid -> weight by cos(el)
    # for solid-angle correction)
    upper = lower = 0.0
    for row in elev_rows:
        w = np.cos(np.radians(row["elevation_deg"]))
        I = 10 ** (row["spl_tower_db"] / 10)
        if row["elevation_deg"] >= 0:
            upper += I * w
        else:
            lower += I * w
    frac_lower = lower / (upper + lower)
    print(f"[{source_case}] power-weighted fraction radiated into lower "
          f"hemisphere (elevation<0): {frac_lower:.2%}")
