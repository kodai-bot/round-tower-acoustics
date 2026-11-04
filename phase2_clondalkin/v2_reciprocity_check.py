"""
Reciprocity check on the v2 interior-A rigid-cavity BEM formulation --
carried forward from ../phase1_generic_tower/review_round_notes.md Task 5 / v1's
check_reciprocity.py, re-run on the new geometry per
rebuild_spec.md SS8.

Acoustic reciprocity: for a fixed domain and boundary condition, the
transfer pressure between two points A and B does not depend on which one
is the source. A cheap, independent check on the interior formulation
(same equation validated in v1's validate_interior_bem.py) that exercises
the ACTUAL v2 tower geometry rather than just the closed-form sphere test.

    (0.5*Id + K) u = u_inc|Gamma
    u(x) = u_inc(x) - (D u)(x)      [x strictly interior]

Points: v1 used (0,0,1.75) and (1.9,0,0.7) -- both OUTSIDE v2's much
smaller r_inner=0.9905m chamber, so new points are used, sized to this
geometry: the source position actually used elsewhere in this round
(chamber centre, window-band height) and a point near window 0, still
strictly interior.

Test frequency: 150Hz (the project's existing SPL-calibration reference
frequency, REFERENCE_FREQ_HZ) -- a representative, physically meaningful
choice, not tied to whichever frequency v2's own production sweep later
turns out to favour (that result isn't known yet when this gate-adjacent
check runs, and reciprocity should hold at any frequency if the
formulation is correct, not just a "best" one).
"""
import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential

import v2_params as P

Q = 0.01731   # arbitrary fixed monopole strength (cancels out in the
              # magnitude ratio check below), same value v1 used
TEST_FREQ = P.REFERENCE_FREQ_HZ

POINT_A = np.array(P.SOURCE_POSITION)         # chamber centre, window-band height
POINT_B = np.array([0.85, 0.0, 1.7])          # near window 0 (r=0.85 < r_inner=0.9905,
                                                # window 0 bearing, within window
                                                # z-band 1.219-2.311), strictly interior

grid = bempp.import_grid(P.MESH_INTERIOR_A_PATH)
space = bempp.function_space(grid, "DP", 0)
identity = bempp.operators.boundary.sparse.identity(space, space, space)

omega = 2 * np.pi * TEST_FREQ
k_real = omega / P.C_SOUND
k = k_real * (1 - 1j * P.LOSS_DELTA)

dlp = bempp.operators.boundary.helmholtz.double_layer(space, space, space, k)
lhs_op = 0.5 * identity + dlp


def u_inc_at(x, source):
    dx = x - source
    r = max(np.linalg.norm(dx), 1e-9)
    return 1j * omega * P.RHO0 * Q / (4 * np.pi * r) * np.exp(1j * k * r)


def solve_and_evaluate(source, receiver):
    # NOTE: numba's nopython JIT (used by bempp's @complex_callable) can't
    # call back into an arbitrary Python function -- the u_inc formula has
    # to be inlined directly here rather than calling u_inc_at(), exactly
    # the same constraint v1's solve_tower_radiation.py already documents
    # for its window-bucketing logic.
    @bempp.complex_callable
    def dirichlet_data(x, n, domain_index, result):
        dx = x - source
        r = max(np.linalg.norm(dx), 1e-9)
        result[0] = 1j * omega * P.RHO0 * Q / (4 * np.pi * r) * np.exp(1j * k * r)

    uinc_fun = bempp.GridFunction(space, fun=dirichlet_data)
    u_sol, info = gmres(lhs_op, uinc_fun, tol=1e-5)
    if info != 0:
        print(f"  [WARNING] gmres info={info}")

    pt = np.array(receiver).reshape(3, 1)
    dlp_pot = helmholtz_potential.double_layer(space, pt, k)
    p_at_receiver = u_inc_at(receiver, source) - dlp_pot.evaluate(u_sol).flatten()[0]
    return p_at_receiver


if __name__ == "__main__":
    print(f"v2 reciprocity check at {TEST_FREQ} Hz, mesh {P.MESH_INTERIOR_A_PATH}")
    print(f"  Point A (source position, chamber centre): {POINT_A}")
    print(f"  Point B (near window 0):                   {POINT_B}")

    p_forward = solve_and_evaluate(POINT_A, POINT_B)
    p_reversed = solve_and_evaluate(POINT_B, POINT_A)

    print(f"\nsource=A, receiver=B: p={p_forward:.6g}  |p|={abs(p_forward):.5f}  "
          f"phase={np.angle(p_forward):.4f}")
    print(f"source=B, receiver=A: p={p_reversed:.6g}  |p|={abs(p_reversed):.5f}  "
          f"phase={np.angle(p_reversed):.4f}")

    rel_err_mag = abs(abs(p_forward) - abs(p_reversed)) / abs(p_forward)
    phase_diff = abs(np.angle(p_forward) - np.angle(p_reversed))
    phase_diff = min(phase_diff, 2 * np.pi - phase_diff)

    print(f"\nrelative magnitude difference: {rel_err_mag:.4%}")
    print(f"phase difference: {phase_diff:.4f} rad ({np.degrees(phase_diff):.2f} deg)")

    TOL_MAG = 0.05
    TOL_PHASE_DEG = 5.0
    passed = rel_err_mag < TOL_MAG and np.degrees(phase_diff) < TOL_PHASE_DEG
    if passed:
        print(f"\nPASS -- reciprocity holds within {TOL_MAG:.0%} magnitude / "
              f"{TOL_PHASE_DEG}deg phase tolerance (same tolerance v1 used).")
    else:
        print(f"\nFAIL -- reciprocity violated beyond tolerance on the v2 geometry.")
        raise SystemExit(1)
