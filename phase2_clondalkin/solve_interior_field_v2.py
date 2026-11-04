"""
v2 interior sound field solve -- Variant A (ceiled chamber), Clondalkin
dimensions. See v2_params.py for geometry and rebuild_spec.md
for the source data.

Formulation is IDENTICAL to v1's solve_interior_field.py, reused
unchanged (validated in validate_interior_bem.py against the closed-form
monopole-in-sphere solution, and re-confirmed on this project's actual
tower topology via reciprocity -- see v2_reciprocity_check.py for the v2
re-run of that check):

    (0.5*Id + K) u = u_inc|Gamma
    u(x) = u_inc(x) - (D u)(x)

Only the geometry, mesh, and source position differ from v1:
  - Mesh: v2_mesh_interior_A.msh -- a genuinely closed r_inner=0.9905m
    canister (real floor AND ceiling), not v1's zero-thickness
    cylinder+cone shell (which -- see v2_check_watertight.py's docstring
    -- turned out to be an open tube missing a floor).
  - Source: (0,0,1.5m), INSIDE the window band (1.219-2.311m) per the
    instruction's SS3, unlike v1's source which sat above the window tops.
  - Windows: same weak-coupling approximation as v1 -- treated as
    closed/rigid for this step (they're a much larger fraction of the
    shell here, ~11% by area vs v1's ~5%, since the interior mesh is
    smaller; still treated as negligible per the instruction's SS3/SS5,
    which does not revise this simplification -- it remains on the open-
    items list, same as it was in v1).

Only ONE source case this round (the standing, centred, "room as a 4-way
speaker" test) -- v1's seated/near-window control case is not repeated
here; v1 already established the qualitative "seated = lopsided" result
once, at the coarser mesh, and this round's changes don't bear on that
control case's own conclusion.
"""
import csv

import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres

import v2_params as P

OUT_CSV = P.DRIVING_CSV_PATH

# ---------------------------------------------------------------------
p_rms_1m = 20e-6 * 10 ** (P.TARGET_SPL_AT_1M_DB / 20)
p_peak_1m = p_rms_1m * np.sqrt(2)
omega_ref = 2 * np.pi * P.REFERENCE_FREQ_HZ
Q = p_peak_1m * 4 * np.pi / (omega_ref * P.RHO0)
print(f"target {P.TARGET_SPL_AT_1M_DB}dB SPL at 1m @ {P.REFERENCE_FREQ_HZ}Hz "
      f"-> constant source volume velocity Q={Q:.5f} m^3/s")

grid = bempp.import_grid(P.MESH_INTERIOR_A_PATH)
space = bempp.function_space(grid, "DP", 0)
identity = bempp.operators.boundary.sparse.identity(space, space, space)

domain_indices = grid.domain_indices
centroids = grid.centroids
volumes = grid.volumes

window_mask = domain_indices == P.WINDOWS_DOMAIN
window_azimuth = np.degrees(np.arctan2(centroids[:, 1], centroids[:, 0])) % 360
window_ref_angles = np.array([i * 360 / P.NUM_WINDOWS for i in range(P.NUM_WINDOWS)])


def nearest_window_index(az_deg):
    diffs = np.abs((az_deg - window_ref_angles + 180) % 360 - 180)
    return int(np.argmin(diffs))


window_index = np.full(grid.number_of_elements, -1, dtype=int)
for i in np.where(window_mask)[0]:
    window_index[i] = nearest_window_index(window_azimuth[i])

for w in range(P.NUM_WINDOWS):
    n_elem = np.sum(window_index == w)
    area = volumes[window_index == w].sum()
    print(f"window {w} (azimuth {window_ref_angles[w]:.0f} deg): "
          f"{n_elem} elements, {area:.3f} m^2")

total_shell_area = volumes.sum()
window_area_total = volumes[window_mask].sum()
print(f"window area fraction of shell: {window_area_total / total_shell_area:.1%} "
      f"(v1's interior chamber was ~5%; this v2 canister is smaller so the "
      f"fraction is higher -- still treated as negligible for this "
      f"weak-coupling step, per the instruction's SS3/SS5)")


def solve_case(f, source_pos):
    """(0.5*Id + K) u = u_inc|Gamma, TOTAL boundary pressure trace u
    directly -- see v1's validate_interior_bem.py / this project's
    reciprocity check for validation of this exact formulation."""
    omega = 2 * np.pi * f
    k_real = omega / P.C_SOUND
    k = k_real * (1 - 1j * P.LOSS_DELTA)

    dlp = bempp.operators.boundary.helmholtz.double_layer(space, space, space, k)
    lhs_op = 0.5 * identity + dlp

    @bempp.complex_callable
    def dirichlet_data(x, n, domain_index, result):
        dx = x - np.array(source_pos)
        r = max(np.linalg.norm(dx), 1e-6)
        result[0] = 1j * omega * P.RHO0 * Q / (4 * np.pi * r) * np.exp(1j * k * r)

    uinc_fun = bempp.GridFunction(space, fun=dirichlet_data)
    u_sol, info = gmres(lhs_op, uinc_fun, tol=1e-5, maxiter=P.GMRES_MAXITER)
    if info != 0:
        print(f"  [WARNING] f={f}Hz source={source_pos}: gmres did not "
              f"cleanly converge within {P.GMRES_MAXITER} iterations "
              f"(info={info})", flush=True)
    return u_sol.coefficients


if __name__ == "__main__":
    import time
    rows = []
    sweep_start = time.time()
    for i, f in enumerate(P.FREQUENCIES):
        print(f"solving interior field: f={f}Hz ... ({i}/{len(P.FREQUENCIES)} done, "
              f"{time.time()-sweep_start:.0f}s elapsed)", flush=True)
        u_boundary = solve_case(f, P.SOURCE_POSITION)
        for w in range(P.NUM_WINDOWS):
            mask = window_index == w
            area = volumes[mask]
            p_window = np.sum(u_boundary[mask] * area) / area.sum()
            rows.append({
                "source_case": "standing",
                "freq_hz": f,
                "window_index": w,
                "azimuth_deg": window_ref_angles[w],
                "pressure_re": p_window.real,
                "pressure_im": p_window.imag,
                "pressure_mag": abs(p_window),
            })

        # checkpoint every frequency -- see solve_tower_radiation_v2.py's
        # matching comment for why (real-time progress visibility + partial
        # results survive a kill)
        with open(OUT_CSV, "w", newline="") as fh:
            for line in P.provenance_header(
                    "solve_interior_field_v2.py",
                    extra={"mesh": P.MESH_INTERIOR_A_PATH, "source_position": P.SOURCE_POSITION}):
                fh.write(line + "\n")
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    print(f"\nwrote {OUT_CSV} ({len(rows)} rows), total {time.time()-sweep_start:.0f}s", flush=True)

    print(f"\n{'freq':>6} {'min_mag':>10} {'max_mag':>10} {'imbalance_dB':>13}")
    for f in P.FREQUENCIES:
        mags = [r["pressure_mag"] for r in rows if r["freq_hz"] == f]
        imbalance_db = 20 * np.log10(max(mags) / max(min(mags), 1e-12))
        print(f"{f:6.0f} {min(mags):10.5f} {max(mags):10.5f} {imbalance_db:13.2f}")
