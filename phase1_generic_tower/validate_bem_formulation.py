"""
Validation test: solve the exterior Helmholtz Neumann (prescribed-velocity)
radiation problem on a sphere using the direct BIE

    (0.5*Id - K) u = -V g

and compare against the exact analytical solution for a pulsating sphere
(a uniformly breathing sphere of radius a with prescribed radial surface
velocity amplitude U0), which is a standard textbook result (Kinsler &
Frey / Fahy):

    p(r) = rho0 * c * U0 * (i*k*a)/(1 + i*k*a) * (a/r) * exp(i*k*(r-a))

This lets us verify the sign convention of the boundary integral equation
against ground truth before trusting it on the actual tower geometry.
"""
import numpy as np
import bempp_cl.api as bempp

c = 343.0
rho0 = 1.225
f = 150.0
k = 2 * np.pi * f / c

a = 1.0       # sphere radius
U0 = 1.0      # uniform prescribed radial velocity amplitude (real)

grid = bempp.shapes.regular_sphere(4)
space = bempp.function_space(grid, "DP", 0)

identity = bempp.operators.boundary.sparse.identity(space, space, space)
dlp = bempp.operators.boundary.helmholtz.double_layer(space, space, space, k)
slp = bempp.operators.boundary.helmholtz.single_layer(space, space, space, k)

lhs = 0.5 * identity - dlp


omega = 2 * np.pi * f
# BIE Neumann data is dp/dn, related to prescribed particle velocity via
# the linearized Euler equation (dp/dn = i*omega*rho0*v_n for e^{-i*omega*t}
# convention): NOT simply U0 itself.
g_scale = 1j * omega * rho0 * U0


@bempp.complex_callable
def neumann_data(x, n, domain_index, result):
    result[0] = g_scale


g_fun = bempp.GridFunction(space, fun=neumann_data)
rhs = -1.0 * (slp * g_fun)

from bempp_cl.api.linalg import gmres
u_sol, info = gmres(lhs, rhs, tol=1e-6)
print("gmres info:", info)

# Evaluate at a handful of exterior points along +x axis
test_r = np.array([1.5, 2.0, 3.0, 5.0, 10.0])
points = np.vstack([test_r, np.zeros_like(test_r), np.zeros_like(test_r)])

from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential
dlp_pot = helmholtz_potential.double_layer(space, points, k)
slp_pot = helmholtz_potential.single_layer(space, points, k)

p_bem = (dlp_pot.evaluate(u_sol) - slp_pot.evaluate(g_fun)).flatten()

p_exact = rho0 * c * U0 * (a / test_r) * (1j * k * a) / (1j * k * a - 1) * np.exp(1j * k * (test_r - a))

print(f"{'r':>6} {'|p_bem|':>12} {'|p_exact|':>12} {'phase_bem':>12} {'phase_exact':>12}")
for r, pb, pe in zip(test_r, p_bem, p_exact):
    print(f"{r:6.2f} {abs(pb):12.5f} {abs(pe):12.5f} {np.angle(pb):12.5f} {np.angle(pe):12.5f}")

rel_err = np.abs(p_bem - p_exact) / np.abs(p_exact)
print("relative errors:", rel_err)
print("max relative error:", rel_err.max())
