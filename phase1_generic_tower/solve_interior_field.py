"""
Interior sound field of the Round Tower upper chamber, driven by an
omnidirectional point-source monopole (standing in for a chanting voice),
evaluated at each of the 4 window locations individually.

This directly tests the "room as a 4-way speaker" hypothesis: does a
centred, unaimed source get redistributed roughly evenly across all 4
windows by the room's own shape (walls + dome), or does it stay
lopsided? It replaces the old flat WINDOW_SPL_DB guess in
solve_tower_radiation.py with 4 individually computed driving values.

Two source cases (see project chat log for the reasoning):
  - "standing": centre of the room, mouth height ~1.75m. This is well
    above the window tops (1.45m) -- there's no direct line-of-sight
    aiming at any window from a natural standing posture, so whatever
    reaches the windows has to get there via wall/dome reflection. This
    is the actual hypothesis under test.
  - "seated": near window 0, roughly window-centre height, as a control
    case reproducing "sitting and singing straight out one window" --
    known from real bell tests to show no effect. If this model shows a
    similarly muted/unremarkable exterior result, that's a sanity check
    that the model isn't producing nonsense.

Formulation: rigid closed cavity (windows treated as closed for this
step -- standard weak-coupling approximation, they're ~5% of the total
shell area), validated in validate_interior_bem.py:

    (0.5*Id + K) u = u_inc|Gamma
    u(x) = u_inc(x) - (D u)(x)

u_inc is the point source's own free-space field. A small uniform loss
factor is folded into the wavenumber (k -> k*(1 - i*delta)) since a
lossless rigid cavity has genuine (not numerical) resonances that would
otherwise blow up -- there's no real absorption measurement for this site
to calibrate against, so this is a placeholder, not a measured value.

Window driving level: rather than guess an SPL at the aperture (as the
earlier version of this pipeline did), the source itself is calibrated to
a target free-field SPL at 1m (a much more physically intuitive
"how loud is the voice" reference point), held at a CONSTANT volume
velocity amplitude across the frequency sweep.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres

MESH_PATH = f"{HERE}/tower_upper_chamber.msh"
OUT_CSV = f"{HERE}/window_driving_conditions.csv"

c = 343.0
rho0 = 1.225

TARGET_SPL_AT_1M_DB = 84.0   # sustained male singing voice at 1m -- revised
                              # down from an earlier 95dB guess (which sat
                              # in "exceptional peaks" territory, >95dB,
                              # not sustainable delivery). 84dB sits between
                              # "amateur singing" (70-80dB) and "trained
                              # classical male singer" (80-90dB) in typical
                              # published SPL-at-1m tables, consistent with
                              # muezzin call-to-prayer measurements (~80-90dB
                              # at 1m). Still a guess -- no measurement
                              # exists for sean-nos/psalm-precentor delivery
                              # specifically -- but a better-grounded one.
REFERENCE_FREQ_HZ = 150.0    # calibration reference frequency
LOSS_DELTA = 0.01            # small placeholder damping (k -> k*(1-i*delta))

WINDOWS_DOMAIN = 1
WALLS_DOMAIN = 2
NUM_WINDOWS = 4
TOWER_RADIUS = 2.2  # must match build_tower_mesh.py

SOURCE_STANDING = (0.0, 0.0, 1.75)          # main test: centred, standing
SOURCE_SEATED = (1.9, 0.0, 0.7)             # control: seated near window 0
                                              # (azimuth 0), window-centre
                                              # height (window spans 0-1.45m)

FREQUENCIES = [50, 75, 100, 125, 143, 150, 160, 165, 170, 175, 180, 185, 190,
               196, 200, 207, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
               310, 320, 330, 340, 349, 360, 375, 390]
# extended to ~390Hz (mesh re-cut to 0.15m elements, see build_tower_mesh.py)
# to cover French's 5th and 6th scale degrees (~294-349Hz), which the
# original 300Hz-ceiling mesh couldn't resolve

# ---------------------------------------------------------------------
p_rms_1m = 20e-6 * 10 ** (TARGET_SPL_AT_1M_DB / 20)
p_peak_1m = p_rms_1m * np.sqrt(2)
omega_ref = 2 * np.pi * REFERENCE_FREQ_HZ
Q = p_peak_1m * 4 * np.pi / (omega_ref * rho0)
print(f"target {TARGET_SPL_AT_1M_DB}dB SPL at 1m @ {REFERENCE_FREQ_HZ}Hz "
      f"-> constant source volume velocity Q={Q:.5f} m^3/s")

grid = bempp.import_grid(MESH_PATH)
space = bempp.function_space(grid, "DP", 0)
identity = bempp.operators.boundary.sparse.identity(space, space, space)

domain_indices = grid.domain_indices
centroids = grid.centroids
volumes = grid.volumes

window_mask = domain_indices == WINDOWS_DOMAIN
window_azimuth = np.degrees(np.arctan2(centroids[:, 1], centroids[:, 0])) % 360
window_ref_angles = np.array([i * 360 / NUM_WINDOWS for i in range(NUM_WINDOWS)])


def nearest_window_index(az_deg):
    diffs = np.abs((az_deg - window_ref_angles + 180) % 360 - 180)
    return int(np.argmin(diffs))


window_index = np.full(grid.number_of_elements, -1, dtype=int)
for i in np.where(window_mask)[0]:
    window_index[i] = nearest_window_index(window_azimuth[i])

for w in range(NUM_WINDOWS):
    n_elem = np.sum(window_index == w)
    area = volumes[window_index == w].sum()
    print(f"window {w} (azimuth {window_ref_angles[w]:.0f} deg): "
          f"{n_elem} elements, {area:.3f} m^2")


def solve_case(f, source_pos):
    """Solves (0.5*Id + K) u = u_inc|Gamma for the TOTAL boundary pressure
    trace u directly (see validate_interior_bem.py) -- u_sol's own
    coefficients already ARE u|Gamma, no further combination needed."""
    omega = 2 * np.pi * f
    k_real = omega / c
    k = k_real * (1 - 1j * LOSS_DELTA)

    dlp = bempp.operators.boundary.helmholtz.double_layer(space, space, space, k)
    lhs_op = 0.5 * identity + dlp

    @bempp.complex_callable
    def dirichlet_data(x, n, domain_index, result):
        dx = x - np.array(source_pos)
        r = max(np.linalg.norm(dx), 1e-6)
        result[0] = 1j * omega * rho0 * Q / (4 * np.pi * r) * np.exp(1j * k * r)

    uinc_fun = bempp.GridFunction(space, fun=dirichlet_data)
    u_sol, info = gmres(lhs_op, uinc_fun, tol=1e-5)
    if info != 0:
        print(f"  [WARNING] f={f}Hz source={source_pos}: gmres info={info}")
    return u_sol.coefficients


rows = []
for source_case, source_pos in [("standing", SOURCE_STANDING)]:
    # NOTE: "seated" control case intentionally skipped for this extended-
    # range (0.15m mesh) run to keep compute time reasonable (~50s/solve
    # here vs a few seconds at the old 0.2m/300Hz mesh). It was already
    # conclusively established as the lopsided/negative case at the
    # coarser mesh; re-add here if it's worth re-confirming at higher
    # resolution later.
    for f in FREQUENCIES:
        print(f"solving interior field: f={f}Hz source={source_case} ...")
        u_boundary = solve_case(f, source_pos)
        for w in range(NUM_WINDOWS):
            mask = window_index == w
            area = volumes[mask]
            p_window = np.sum(u_boundary[mask] * area) / area.sum()
            rows.append({
                "source_case": source_case,
                "freq_hz": f,
                "window_index": w,
                "azimuth_deg": window_ref_angles[w],
                "pressure_re": p_window.real,
                "pressure_im": p_window.imag,
                "pressure_mag": abs(p_window),
            })

with open(OUT_CSV, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
print(f"\nwrote {OUT_CSV} ({len(rows)} rows)")

# quick imbalance summary
print(f"\n{'case':>10} {'freq':>6} {'min_mag':>10} {'max_mag':>10} {'imbalance_dB':>13}")
for source_case in ("standing",):
    for f in FREQUENCIES:
        mags = [r["pressure_mag"] for r in rows if r["source_case"] == source_case and r["freq_hz"] == f]
        imbalance_db = 20 * np.log10(max(mags) / max(min(mags), 1e-12))
        print(f"{source_case:>10} {f:6.0f} {min(mags):10.5f} {max(mags):10.5f} {imbalance_db:13.2f}")
