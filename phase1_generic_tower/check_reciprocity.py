"""
Reciprocity check on the interior rigid-cavity BEM formulation (Task 5,
review_round_notes.md).

Acoustic reciprocity (Rayleigh reciprocity theorem): for a fixed domain and
boundary condition, the transfer pressure between two points A and B does
not depend on which one is the source -- G(A,B) = G(B,A). This is a cheap,
independent check on the interior formulation validated in
validate_interior_bem.py: if source/receiver swap doesn't reproduce the
same transfer pressure, that indicates a real bug (e.g. a sign/jump-
relation error) that the closed-form sphere validation might not catch,
since this checks a genuinely different pair of points on the actual tower
geometry rather than a symmetric test case.

Uses the same validated equation as solve_interior_field.py:
    (0.5*Id + K) u = u_inc|Gamma
    u(x) = u_inc(x) - (D u)(x)      [x strictly interior]

Points used: the same standing-source and seated-source positions already
used elsewhere in this project (SOURCE_STANDING, SOURCE_SEATED in
solve_interior_field.py), reused here purely as two convenient, already-
meaningful interior points -- not because either represents a real singer
in this test.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential

MESH_PATH = f"{HERE}/tower_upper_chamber.msh"

c = 343.0
rho0 = 1.225
Q = 0.01731  # arbitrary fixed monopole strength (cancels out in the ratio
              # check below, but keep it physically reasonable)
LOSS_DELTA = 0.01
TEST_FREQ = 185.0  # the model's best-projecting frequency

POINT_A = np.array([0.0, 0.0, 1.75])   # "standing" position (chamber centre)
POINT_B = np.array([1.9, 0.0, 0.7])    # "seated" position (near window 0)

grid = bempp.import_grid(MESH_PATH)
space = bempp.function_space(grid, "DP", 0)
identity = bempp.operators.boundary.sparse.identity(space, space, space)

omega = 2 * np.pi * TEST_FREQ
k_real = omega / c
k = k_real * (1 - 1j * LOSS_DELTA)

dlp = bempp.operators.boundary.helmholtz.double_layer(space, space, space, k)
lhs_op = 0.5 * identity + dlp


def u_inc(x, source):
    dx = np.atleast_2d(x) - source
    r = np.linalg.norm(dx, axis=-1)
    r = np.maximum(r, 1e-9)
    return 1j * omega * rho0 * Q / (4 * np.pi * r) * np.exp(1j * k * r)


def solve_and_evaluate(source, receiver):
    @bempp.complex_callable
    def dirichlet_data(x, n, domain_index, result):
        dx = x - source
        r = max(np.linalg.norm(dx), 1e-9)
        result[0] = 1j * omega * rho0 * Q / (4 * np.pi * r) * np.exp(1j * k * r)

    uinc_fun = bempp.GridFunction(space, fun=dirichlet_data)
    u_sol, info = gmres(lhs_op, uinc_fun, tol=1e-5)
    if info != 0:
        print(f"  [WARNING] gmres info={info}")

    pt = np.array(receiver).reshape(3, 1)
    dlp_pot = helmholtz_potential.double_layer(space, pt, k)
    p_at_receiver = u_inc(receiver, source)[0] - dlp_pot.evaluate(u_sol).flatten()[0]
    return p_at_receiver


print(f"Reciprocity check at {TEST_FREQ} Hz")
print(f"  Point A (chamber centre): {POINT_A}")
print(f"  Point B (near window 0):  {POINT_B}")

p_forward = solve_and_evaluate(POINT_A, POINT_B)   # source A, receive B
p_reversed = solve_and_evaluate(POINT_B, POINT_A)  # source B, receive A

print(f"\nsource=A, receiver=B: p = {p_forward:.6g}  |p|={abs(p_forward):.5f}  phase={np.angle(p_forward):.4f}")
print(f"source=B, receiver=A: p = {p_reversed:.6g}  |p|={abs(p_reversed):.5f}  phase={np.angle(p_reversed):.4f}")

rel_err_mag = abs(abs(p_forward) - abs(p_reversed)) / abs(p_forward)
phase_diff = abs(np.angle(p_forward) - np.angle(p_reversed))
phase_diff = min(phase_diff, 2 * np.pi - phase_diff)

print(f"\nrelative magnitude difference: {rel_err_mag:.4%}")
print(f"phase difference: {phase_diff:.4f} rad ({np.degrees(phase_diff):.2f} deg)")

TOL_MAG = 0.05
TOL_PHASE_DEG = 5.0
if rel_err_mag < TOL_MAG and np.degrees(phase_diff) < TOL_PHASE_DEG:
    print(f"\nPASS -- reciprocity holds within {TOL_MAG:.0%} magnitude / "
          f"{TOL_PHASE_DEG}deg phase tolerance. No evidence of a sign/jump-"
          f"relation bug in the interior formulation.")
else:
    print(f"\nFAIL -- reciprocity violated beyond tolerance. This indicates "
          f"a real problem with the interior BEM formulation that the "
          f"closed-form sphere validation did not catch.")
