"""
Cross-check on the 250Hz trough/directional-null result (Task 6,
review_round_notes.md): is it a genuine acoustic feature, or an
artifact of the direct (non-CFIE) method's known weakness near "irregular
frequencies" (spurious non-uniqueness at frequencies coinciding with
interior Dirichlet eigenvalues of the enclosed volume)?

METHOD SUBSTITUTION, stated plainly: the instructions suggested a true
CFIE/Burton-Miller combined formulation. A full Burton-Miller regularization
of the DIRECT method (unknown = Dirichlet trace given known Neumann data)
needs the hypersingular operator, which this project has deliberately
avoided elsewhere due to real implementation risk (see project history --
an earlier, simpler interior formulation attempt was subtly wrong before
being caught by validation). Building an unvalidated hypersingular-operator
formulation specifically to check for bugs is a bad trade: a mistake there
could produce a false "confirmation" that's worse than not checking at all.

Instead this uses an INDIRECT single-layer formulation of the SAME physical
problem -- a different integral equation with a DIFFERENT (generally
non-overlapping) set of irregular frequencies than the direct method,
using only operators (V, K') already used and validated elsewhere in this
project. If the two independently-formulated methods, which fail at
different frequencies when they fail at all, agree at 250Hz, that's strong
evidence the result is genuine physics, not a direct-method artifact.

Indirect single-layer ansatz: p(x) = S[phi](x), x exterior, unknown density
phi solving (using the exterior jump relation for the single-layer
potential's normal derivative):
    (-0.5*Id + K') phi = g
Validated below against the same closed-form pulsating-sphere solution
used in validate_bem_formulation.py, BEFORE trusting it on the tower.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential

c = 343.0
rho0 = 1.225

# ---------------------------------------------------------------------
# Step 1: validate the indirect single-layer formulation against the
# pulsating sphere (same test as validate_bem_formulation.py)
# ---------------------------------------------------------------------
print("=== Step 1: validating indirect single-layer formulation ===")
a = 1.0
U0 = 1.0
f_val = 150.0
omega_val = 2 * np.pi * f_val
k_val = omega_val / c

grid_s = bempp.shapes.regular_sphere(4)
space_s = bempp.function_space(grid_s, "DP", 0)
identity_s = bempp.operators.boundary.sparse.identity(space_s, space_s, space_s)
adlp_s = bempp.operators.boundary.helmholtz.adjoint_double_layer(space_s, space_s, space_s, k_val)

g_scale_val = 1j * omega_val * rho0 * U0


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

p_exact = rho0 * c * U0 * (a / test_r) * (1j * k_val * a) / (1j * k_val * a - 1) * np.exp(1j * k_val * (test_r - a))
rel_err = np.abs(p_bem - p_exact) / np.abs(p_exact)
print(f"max relative error vs closed form: {rel_err.max():.4%}")

if rel_err.max() > 0.02:
    raise SystemExit("VALIDATION FAILED -- indirect single-layer formulation "
                      "does not match closed form. Do not trust the 250Hz "
                      "cross-check below.")
print("VALIDATION PASSED -- proceeding to the tower cross-check.\n")

# ---------------------------------------------------------------------
# Step 2: apply to the actual tower mesh at 250Hz, compare against the
# already-computed direct-method result at the same frequency
# ---------------------------------------------------------------------
import csv

MESH_PATH = f"{HERE}/tower_upper_chamber.msh"
DRIVING_CSV = f"{HERE}/window_driving_conditions.csv"
GROUND_CSV = f"{HERE}/ground_projection_standing.csv"

WINDOWS_DOMAIN = 1
NUM_WINDOWS = 4
WINDOW_REF_ANGLES = np.array([i * 360 / NUM_WINDOWS for i in range(NUM_WINDOWS)])
TEST_FREQ = 250.0
TOWER_ELEVATION = 30.0
LISTENER_HEIGHT = 1.6
TEST_AZIMUTHS = [0, 45]  # on-axis (window direction) vs diagonal
TEST_DISTANCE = 100.0

print("=== Step 2: tower cross-check at 250Hz ===")
grid = bempp.import_grid(MESH_PATH)
space = bempp.function_space(grid, "DP", 0)
identity = bempp.operators.boundary.sparse.identity(space, space, space)

driving = {}
with open(DRIVING_CSV) as fh:
    for row in csv.DictReader(fh):
        key = (row["source_case"], int(float(row["freq_hz"])))
        driving.setdefault(key, [None] * NUM_WINDOWS)
        driving[key][int(row["window_index"])] = complex(float(row["pressure_re"]), float(row["pressure_im"]))

p_windows = np.array(driving[("standing", int(TEST_FREQ))], dtype=np.complex128)
ref_angles = WINDOW_REF_ANGLES.copy()

omega = 2 * np.pi * TEST_FREQ
k = omega / c

adlp = bempp.operators.boundary.helmholtz.adjoint_double_layer(space, space, space, k)
lhs = -0.5 * identity + adlp


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
phi, info = gmres(lhs, g_fun, tol=1e-5)
print(f"gmres info: {info}")
if info != 0:
    print("  [WARNING] gmres did not converge cleanly")

pts = []
for az_deg in TEST_AZIMUTHS:
    az = np.radians(az_deg)
    x = TEST_DISTANCE * np.cos(az)
    y = TEST_DISTANCE * np.sin(az)
    z = -TOWER_ELEVATION + LISTENER_HEIGHT
    pts.append((x, y, z))
pts = np.array(pts).T

slp_pot = helmholtz_potential.single_layer(space, pts, k)
p_indirect = slp_pot.evaluate(phi).flatten()


def spl_db(p_complex):
    p_rms = np.abs(p_complex) / np.sqrt(2)
    return 20 * np.log10(np.maximum(p_rms, 1e-12) / 20e-6)


# free-field baseline (same as solve_tower_radiation.py)
window_elems = grid.domain_indices == WINDOWS_DOMAIN
elem_volumes = grid.volumes
window_area_each = float(elem_volumes[window_elems].sum()) / NUM_WINDOWS
Q = sum(abs(p) / (rho0 * c) * window_area_each for p in p_windows)
r_slant = np.sqrt(pts[0] ** 2 + pts[1] ** 2 + pts[2] ** 2)
p_free = 1j * omega * rho0 * Q / (4 * np.pi * r_slant) * np.exp(1j * k * r_slant)
gain_indirect = spl_db(p_indirect) - spl_db(p_free)

# pull the existing direct-method result at the same points for comparison
direct_rows = list(csv.DictReader(open(GROUND_CSV)))
for r in direct_rows:
    for kk in ("freq_hz", "distance_m", "azimuth_deg", "gain_db"):
        r[kk] = float(r[kk])

print(f"\n{'azimuth':>8} {'nearest_avail':>14} {'direct gain_db':>15} {'indirect(SLP) gain_db':>22} {'diff':>8}")
for az_deg, g_ind in zip(TEST_AZIMUTHS, gain_indirect):
    candidates = [r for r in direct_rows
                  if r["freq_hz"] == TEST_FREQ and r["distance_m"] == TEST_DISTANCE]
    nearest = min(candidates, key=lambda r: abs((r["azimuth_deg"] - az_deg + 180) % 360 - 180))
    g_dir = nearest["gain_db"]
    diff = g_ind - g_dir
    print(f"{az_deg:8.0f} {nearest['azimuth_deg']:14.0f} {g_dir:15.2f} {g_ind:22.2f} {diff:8.2f}")

print("\nIf both methods show the same on-axis/diagonal pattern at 250Hz "
      "(on-axis much worse than diagonal), that's evidence the sharp "
      "behavior is genuine physics, not a direct-method-specific "
      "irregular-frequency artifact -- two differently-formulated methods "
      "with different failure modes agreeing is a strong correctness signal.")
