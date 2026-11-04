"""
Direct vs. indirect (CFIE-style) cross-check for v2 -- carried forward
from ../phase1_generic_tower/review_round_notes.md Task 6 / v1's solve_cfie_250hz.py,
per rebuild_spec.md SS8. Re-run at whichever frequency
v2's OWN production sweep shows its sharpest azimuthal feature (v1 used a
hardcoded 250Hz, chosen because that's where v1's headline curve showed
its sharpest trough -- v2's geometry is different enough that the
sharpest feature may sit elsewhere, so this script finds it automatically
from v2_ground_projection_standing.csv rather than assuming).

METHOD SUBSTITUTION, same as v1 (stated plainly there, repeated here):
a full Burton-Miller CFIE needs the hypersingular operator, deliberately
avoided project-wide due to real implementation risk. This instead uses
an INDIRECT single-layer formulation of the same exterior problem -- a
different integral equation with different irregular frequencies than
the direct method, built only from already-validated operators:

    (-0.5*Id + K') phi = g,   p(x) = S[phi](x)

Validated against the closed-form pulsating-sphere solution before being
trusted on the tower (same test as v1's script and
validate_bem_formulation.py).

REQUIRES v2_ground_projection_standing.csv to already exist (produced by
solve_tower_radiation_v2.py's production sweep) -- this script reads it
to pick its target frequency and to pull the direct-method comparison
values, it does not re-run the full sweep itself.
"""
import csv

import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential

import v2_params as P

TEST_DISTANCE = 100.0
TEST_AZIMUTHS = [0, 45]   # on-axis (window bearing) vs diagonal
GROUND_CSV = f"{P.MODEL_DIR}/v2_ground_projection_standing.csv"


def _read_csv_skip_provenance(path):
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    rows = list(csv.DictReader(lines))
    for r in rows:
        for kk in ("freq_hz", "distance_m", "azimuth_deg", "gain_db"):
            r[kk] = float(r[kk])
    return rows


def pick_sharpest_frequency(rows, distance):
    """Frequency (at the given evaluation distance) with the largest
    azimuth-to-azimuth gain_db swing -- i.e. the sharpest directional
    feature, the case most worth stress-testing against a second,
    independently-formulated method."""
    by_freq = {}
    for r in rows:
        if r["distance_m"] != distance:
            continue
        by_freq.setdefault(r["freq_hz"], []).append(r["gain_db"])
    swings = {f: (max(v) - min(v)) for f, v in by_freq.items()}
    best_f = max(swings, key=swings.get)
    print(f"sharpest azimuthal feature: f={best_f:.0f}Hz, swing={swings[best_f]:.2f}dB "
          f"(at {distance:.0f}m)")
    return best_f


# ---------------------------------------------------------------------
# Step 1: validate the indirect single-layer formulation against the
# pulsating sphere
# ---------------------------------------------------------------------
def validate_indirect_formulation():
    print("=== Step 1: validating indirect single-layer formulation ===")
    a = 1.0
    U0 = 1.0
    f_val = 150.0
    omega_val = 2 * np.pi * f_val
    k_val = omega_val / P.C_SOUND

    grid_s = bempp.shapes.regular_sphere(4)
    space_s = bempp.function_space(grid_s, "DP", 0)
    identity_s = bempp.operators.boundary.sparse.identity(space_s, space_s, space_s)
    adlp_s = bempp.operators.boundary.helmholtz.adjoint_double_layer(space_s, space_s, space_s, k_val)

    g_scale_val = 1j * omega_val * P.RHO0 * U0

    @bempp.complex_callable
    def neumann_data_s(x, n, domain_index, result):
        result[0] = g_scale_val

    g_fun_s = bempp.GridFunction(space_s, fun=neumann_data_s)
    lhs_s = -0.5 * identity_s + adlp_s
    phi_s, info_s = gmres(lhs_s, g_fun_s, tol=1e-6)
    print(f"gmres info: {info_s}")

    test_r = np.array([1.5, 2.0, 3.0, 5.0, 10.0])
    points_s = np.vstack([test_r, np.zeros_like(test_r), np.zeros_like(test_r)])
    slp_pot_s = helmholtz_potential.single_layer(space_s, points_s, k_val)
    p_bem = slp_pot_s.evaluate(phi_s).flatten()

    p_exact = P.RHO0 * P.C_SOUND * U0 * (a / test_r) * (1j * k_val * a) / (1j * k_val * a - 1) * np.exp(1j * k_val * (test_r - a))
    rel_err = np.abs(p_bem - p_exact) / np.abs(p_exact)
    print(f"max relative error vs closed form: {rel_err.max():.4%}")

    if rel_err.max() > 0.02:
        raise SystemExit("VALIDATION FAILED -- indirect single-layer formulation "
                          "does not match closed form. Do not trust the cross-check below.")
    print("VALIDATION PASSED -- proceeding to the tower cross-check.\n")


# ---------------------------------------------------------------------
# Step 2: apply to the actual v2 exterior mesh at the chosen frequency
# ---------------------------------------------------------------------
def run_crosscheck(test_freq, direct_rows):
    print(f"=== Step 2: tower cross-check at {test_freq:.0f}Hz ===")
    grid = bempp.import_grid(P.MESH_EXTERIOR_PATH)
    space = bempp.function_space(grid, "DP", 0)
    identity = bempp.operators.boundary.sparse.identity(space, space, space)

    driving = {}
    with open(P.DRIVING_CSV_PATH) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    for row in csv.DictReader(lines):
        key = (row["source_case"], int(float(row["freq_hz"])))
        driving.setdefault(key, [None] * P.NUM_WINDOWS)
        driving[key][int(row["window_index"])] = complex(float(row["pressure_re"]), float(row["pressure_im"]))

    p_windows = np.array(driving[("standing", int(test_freq))], dtype=np.complex128)
    ref_angles = np.array([i * 360 / P.NUM_WINDOWS for i in range(P.NUM_WINDOWS)])

    omega = 2 * np.pi * test_freq
    k = omega / P.C_SOUND

    adlp = bempp.operators.boundary.helmholtz.adjoint_double_layer(space, space, space, k)
    lhs = -0.5 * identity + adlp

    @bempp.complex_callable
    def neumann_data(x, n, domain_index, result):
        if domain_index != P.WINDOWS_DOMAIN:
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
    phi, info = gmres(lhs, g_fun, tol=1e-5)
    print(f"gmres info: {info}")
    if info != 0:
        print("  [WARNING] gmres did not converge cleanly")

    pts = []
    for az_deg in TEST_AZIMUTHS:
        az = np.radians(az_deg)
        x = TEST_DISTANCE * np.cos(az)
        y = TEST_DISTANCE * np.sin(az)
        z = -P.TOWER_ELEVATION + P.LISTENER_HEIGHT
        pts.append((x, y, z))
    pts = np.array(pts).T

    slp_pot = helmholtz_potential.single_layer(space, pts, k)
    p_indirect = slp_pot.evaluate(phi).flatten()

    window_elems = grid.domain_indices == P.WINDOWS_DOMAIN
    elem_volumes = grid.volumes
    window_area_each = float(elem_volumes[window_elems].sum()) / P.NUM_WINDOWS
    Q = sum(abs(p) / (P.RHO0 * P.C_SOUND) * window_area_each for p in p_windows)
    r_slant = np.sqrt(pts[0] ** 2 + pts[1] ** 2 + pts[2] ** 2)
    p_free = 1j * omega * P.RHO0 * Q / (4 * np.pi * r_slant) * np.exp(1j * k * r_slant)

    def spl_db(p_complex):
        p_rms = np.abs(p_complex) / np.sqrt(2)
        return 20 * np.log10(np.maximum(p_rms, 1e-12) / 20e-6)

    gain_indirect = spl_db(p_indirect) - spl_db(p_free)

    print(f"\n{'azimuth':>8} {'nearest_avail':>14} {'direct gain_db':>15} "
          f"{'indirect(SLP) gain_db':>22} {'diff':>8}")
    for az_deg, g_ind in zip(TEST_AZIMUTHS, gain_indirect):
        candidates = [r for r in direct_rows
                      if r["freq_hz"] == test_freq and r["distance_m"] == TEST_DISTANCE]
        nearest = min(candidates, key=lambda r: abs((r["azimuth_deg"] - az_deg + 180) % 360 - 180))
        g_dir = nearest["gain_db"]
        diff = g_ind - g_dir
        print(f"{az_deg:8.0f} {nearest['azimuth_deg']:14.0f} {g_dir:15.2f} {g_ind:22.2f} {diff:8.2f}")

    print("\nIf both methods show the same on-axis/diagonal pattern at this "
          "frequency, that's evidence the sharp directional behavior is "
          "genuine physics, not a direct-method-specific irregular-frequency "
          "artifact -- two differently-formulated methods with different "
          "failure modes agreeing is a strong correctness signal. Per v1's "
          "own experience at 250Hz, expect qualitative agreement but not "
          "necessarily tight quantitative agreement in a numerically "
          "delicate zone -- report both, don't force a verdict.")


if __name__ == "__main__":
    direct_rows = _read_csv_skip_provenance(GROUND_CSV)
    test_freq = pick_sharpest_frequency(direct_rows, TEST_DISTANCE)
    validate_indirect_formulation()
    run_crosscheck(test_freq, direct_rows)
