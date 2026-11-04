"""
v2 exterior acoustic radiation solve -- Variant A (ceiled chamber),
Clondalkin dimensions, nominal wall thickness (0.745m). See v2_params.py
and rebuild_spec.md.

Formulation is IDENTICAL to v1's solve_tower_radiation.py, reused
unchanged (validated against the exact pulsating-sphere solution in
validate_bem_formulation.py):

    (0.5*Id - K) p = -V g
    p(x) = D[p](x) - S[g](x)   for exterior x

What's different from v1:
  - Mesh: v2_mesh_exterior.msh -- a real ~0.745m-thick wall with the 4
    windows modelled as rectangular DUCTS (occ.cut through the full
    thickness), not v1's flush zero-thickness patches.
  - Driving boundary condition: v1 drove the window patch directly (it
    WAS the outer surface, since the shell had no thickness). v2 drives
    the piston at the duct-BACK face -- the flat patch sitting at
    r=r_inner, the INNER mouth of each duct (physical group "windows",
    domain 1; see build_tower_mesh_v2.py's classify_exterior_face()).
    This means the BEM solve itself has to resolve sound propagating the
    full 0.745m through the rigid duct sleeve before radiating from the
    duct's true outer mouth -- this is the entire point of modelling a
    real duct instead of a flush patch; driving at the outer mouth
    instead would make the duct geometry acoustically inert and silently
    collapse v2 back to v1's flush-patch approximation. This bridge is
    still a WEAKLY COUPLED two-step approximation overall (same caveat
    v1 documented): the interior solve (solve_interior_field_v2.py)
    computed this driving pressure treating the window/duct-back
    footprint as closed/rigid, ignoring the ~11% of the interior shell
    area it occupies in this smaller v2 chamber.
  - Evaluation grid: 36 azimuths from the start (not widened mid-round,
    as v1 had to do after finding its original 8-point grid aliased with
    the tower's 4-fold symmetry), plus 3 extra near ground distances
    (15/25/30m, the real Clonmacnois/Kells/Glendalough tower-to-church
    distances) and an extended frequency sweep to 500Hz.
"""
import csv
import time

import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential

import v2_params as P

OUT_DIR = P.MODEL_DIR
WINDOW_REF_ANGLES = np.array([i * 360 / P.NUM_WINDOWS for i in range(P.NUM_WINDOWS)])

# ---------------------------------------------------------------------
grid = bempp.import_grid(P.MESH_EXTERIOR_PATH)
space = bempp.function_space(grid, "DP", 0)
identity = bempp.operators.boundary.sparse.identity(space, space, space)

window_elems = grid.domain_indices == P.WINDOWS_DOMAIN
wall_elems = grid.domain_indices == P.WALLS_DOMAIN
print(f"mesh: {grid.number_of_elements} elements "
      f"({window_elems.sum()} window/duct-back, {wall_elems.sum()} wall/cone/sleeve/floor)")
elem_volumes = grid.volumes
window_area_total = float(elem_volumes[window_elems].sum())
window_area_each = window_area_total / P.NUM_WINDOWS


def _read_driving_csv(path):
    """Skips the leading '# key: value' provenance header lines written
    by v2_params.provenance_header() before handing off to DictReader."""
    driving = {}
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    for row in csv.DictReader(lines):
        key = (row["source_case"], int(float(row["freq_hz"])))
        driving.setdefault(key, [None] * P.NUM_WINDOWS)
        w = int(row["window_index"])
        driving[key][w] = complex(float(row["pressure_re"]), float(row["pressure_im"]))
    return driving


driving = _read_driving_csv(P.DRIVING_CSV_PATH)


def solve_at_frequency(f, source_case="standing"):
    omega = 2 * np.pi * f
    k = omega / P.C_SOUND

    dlp = bempp.operators.boundary.helmholtz.double_layer(space, space, space, k)
    slp = bempp.operators.boundary.helmholtz.single_layer(space, space, space, k)
    lhs = 0.5 * identity - dlp

    p_windows = np.array(driving[(source_case, f)], dtype=np.complex128)
    ref_angles = WINDOW_REF_ANGLES.copy()

    @bempp.complex_callable
    def neumann_data(x, n, domain_index, result):
        if domain_index != P.WINDOWS_DOMAIN:
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

    p_sol, info = gmres(lhs, rhs, tol=1e-5, maxiter=P.GMRES_MAXITER)
    if info != 0:
        print(f"  [WARNING] f={f}Hz source={source_case}: gmres did not "
              f"cleanly converge within {P.GMRES_MAXITER} iterations "
              f"(info={info}) -- possible irregular frequency or "
              f"resolution-margin issue (see v2_params.MESH_SIZE), "
              f"treat result with suspicion", flush=True)

    Q = sum(abs(p) / (P.RHO0 * P.C_SOUND) * window_area_each for p in p_windows)
    return k, omega, p_sol, g_fun, Q


def evaluate_ground_points(k, p_sol, g_fun, Q):
    pts, meta = [], []
    for d in P.GROUND_DISTANCES:
        for az_deg in P.AZIMUTHS_DEG:
            az = np.radians(az_deg)
            x = d * np.cos(az)
            y = d * np.sin(az)
            z = -P.TOWER_ELEVATION + P.LISTENER_HEIGHT
            pts.append((x, y, z))
            meta.append((d, az_deg))
    pts = np.array(pts).T

    dlp_pot = helmholtz_potential.double_layer(space, pts, k)
    slp_pot = helmholtz_potential.single_layer(space, pts, k)
    p_tower = (dlp_pot.evaluate(p_sol) - slp_pot.evaluate(g_fun)).flatten()

    r_slant = np.sqrt(pts[0] ** 2 + pts[1] ** 2 + pts[2] ** 2)
    omega = k * P.C_SOUND
    p_free = 1j * omega * P.RHO0 * Q / (4 * np.pi * r_slant) * np.exp(1j * k * r_slant)
    return meta, p_tower, p_free


def evaluate_elevation_sweep(k, p_sol, g_fun, Q):
    pts, meta = [], []
    for az_deg in P.ELEVATION_AZIMUTHS_DEG:
        for el_deg in P.ELEVATION_ANGLES_DEG:
            az = np.radians(az_deg)
            el = np.radians(el_deg)
            horiz = P.ELEVATION_RADIUS * np.cos(el)
            x = horiz * np.cos(az)
            y = horiz * np.sin(az)
            z = P.ELEVATION_RADIUS * np.sin(el)
            pts.append((x, y, z))
            meta.append((az_deg, el_deg))
    pts = np.array(pts).T

    dlp_pot = helmholtz_potential.double_layer(space, pts, k)
    slp_pot = helmholtz_potential.single_layer(space, pts, k)
    p_tower = (dlp_pot.evaluate(p_sol) - slp_pot.evaluate(g_fun)).flatten()

    omega = k * P.C_SOUND
    p_free = 1j * omega * P.RHO0 * Q / (4 * np.pi * P.ELEVATION_RADIUS) * np.exp(1j * k * P.ELEVATION_RADIUS)
    p_free_arr = np.full_like(p_tower, p_free)
    return meta, p_tower, p_free_arr


def evaluate_pilgrim_rings(k, p_sol, g_fun, Q):
    pts, meta = [], []
    for ring_name, R in P.RING_RADII.items():
        thetas = np.linspace(0, 360, P.N_POINTS_AROUND_RING, endpoint=False)
        site_x = R * np.cos(np.radians(thetas))
        site_y = R * np.sin(np.radians(thetas))
        local_x = site_x + P.TOWER_OFFSET
        local_y = site_y
        local_z = np.full_like(local_x, -P.TOWER_ELEVATION + P.LISTENER_HEIGHT)
        dist_from_tower = np.sqrt(local_x ** 2 + local_y ** 2)
        for th, d, lx, ly, lz in zip(thetas, dist_from_tower, local_x, local_y, local_z):
            pts.append((lx, ly, lz))
            meta.append((ring_name, th, d))
    pts = np.array(pts).T

    dlp_pot = helmholtz_potential.double_layer(space, pts, k)
    slp_pot = helmholtz_potential.single_layer(space, pts, k)
    p_tower = (dlp_pot.evaluate(p_sol) - slp_pot.evaluate(g_fun)).flatten()

    r_slant = np.sqrt(pts[0] ** 2 + pts[1] ** 2 + pts[2] ** 2)
    omega = k * P.C_SOUND
    p_free = 1j * omega * P.RHO0 * Q / (4 * np.pi * r_slant) * np.exp(1j * k * r_slant)
    return meta, p_tower, p_free


def spl_db(p_complex):
    p_rms = np.abs(p_complex) / np.sqrt(2)
    return 20 * np.log10(np.maximum(p_rms, 1e-12) / 20e-6)


def _write(rows, path, extra_prov=None):
    with open(path, "w", newline="") as fh:
        for line in P.provenance_header("solve_tower_radiation_v2.py", extra=extra_prov):
            fh.write(line + "\n")
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    source_case = "standing"
    ground_rows, elev_rows, ring_rows = [], [], []
    sweep_start = time.time()
    n_done = 0
    for f in P.FREQUENCIES:
        if (source_case, f) not in driving:
            print(f"skipping f={f}: no driving data", flush=True)
            continue
        f_start = time.time()
        print(f"solving f={f} Hz ... ({n_done}/{len(P.FREQUENCIES)} done, "
              f"{time.time()-sweep_start:.0f}s elapsed)", flush=True)
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

        n_done += 1
        print(f"  f={f}Hz done in {time.time()-f_start:.1f}s", flush=True)

        # Checkpoint after every frequency, not just at the end -- a full
        # sweep of this size has already once run far longer than
        # estimated (see v2_params.GMRES_MAXITER's comment); overwriting
        # these each iteration means a kill/crash mid-sweep still leaves
        # real partial results and progress is observable in real time
        # from the file's row count, not just inferred from process CPU
        # usage with no visibility into how many frequencies are done.
        _write(ground_rows, f"{OUT_DIR}/v2_ground_projection_standing.csv")
        _write(elev_rows, f"{OUT_DIR}/v2_elevation_sweep_standing.csv")
        _write(ring_rows, f"{OUT_DIR}/v2_pilgrim_ring_walk_standing.csv")

    print(f"wrote all 3 CSVs ({len(ground_rows)}/{len(elev_rows)}/{len(ring_rows)} rows), "
          f"total sweep time {time.time()-sweep_start:.0f}s", flush=True)

    upper = lower = 0.0
    for row in elev_rows:
        w = np.cos(np.radians(row["elevation_deg"]))
        I = 10 ** (row["spl_tower_db"] / 10)
        if row["elevation_deg"] >= 0:
            upper += I * w
        else:
            lower += I * w
    frac_lower = lower / (upper + lower)
    print(f"power-weighted fraction radiated into lower hemisphere "
          f"(elevation<0): {frac_lower:.2%}")
