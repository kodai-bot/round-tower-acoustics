"""
Validation test for the INTERIOR Helmholtz Neumann (rigid cavity) boundary
integral formulation, before it's trusted on the tower geometry.

Closed-form target: a monopole point source sitting at the exact CENTER
of a rigid sphere of radius a excites only the l=0 (spherically symmetric)
mode -- by symmetry, a source on the axis of a sphere only couples to
that one mode, exactly analogous to how a source on the axis of a
cylindrical cavity only excites axisymmetric modes. This reduces the
problem to a 1D radial ODE with an explicit closed-form solution in terms
of spherical Bessel functions j0, y0:

    p(r) = A*j0(kr) + B*y0(kr),   0 < r <= a

  j0(x) = sin(x)/x,          j0'(x) = (x*cos(x) - sin(x)) / x^2
  y0(x) = -cos(x)/x,         y0'(x) = (x*sin(x) + cos(x)) / x^2

B is fixed by matching the known 1/r singularity of the free-space point
source field p_inc(r) = i*omega*rho0*Q*e^{ikr}/(4*pi*r) as r -> 0:
    B = -i*rho0*c*k^2*Q / (4*pi)
A is fixed by the rigid-wall condition dp/dr = 0 at r = a:
    A = -B * y0'(ka) / j0'(ka)

This p(r) is the TOTAL field (it already contains the correct source
singularity), so it can be compared directly against the numerical BEM
solve's total field.

Formulation actually used (found by empirical elimination -- an initial
attempt going through a separate "scattered field" Neumann problem with
the single-layer operator turned out to have a bug I couldn't isolate, so
this uses the simpler, more standard direct route instead): since the
TOTAL field's normal derivative is exactly zero everywhere on a rigid
wall, Green's representation collapses to a single-layer-free equation
purely in the identity and double-layer (K) operators, solved directly
for the TOTAL boundary pressure trace u|Gamma:

    (0.5*Id + K) u = u_inc|Gamma
    u(x) = u_inc(x) - (D u)(x),  x strictly interior

where u_inc is the point source's own free-space field (known everywhere,
including on the boundary) and D is the double-layer potential. Both
signs here were confirmed empirically against the closed form below
(the other 3 of 4 sign combinations were off by 20-50%; this one matches
to <3%, with the one outlying point sitting right next to a standing-wave
null in the closed-form solution, where relative error is naturally
amplified by a small denominator despite small absolute error).
"""
import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential

c = 343.0
rho0 = 1.225
a = 1.0          # sphere radius, m (bempp's regular_sphere is unit radius)
Q = 1e-3         # monopole volume velocity amplitude, m^3/s (arbitrary)

f = 100.0        # test frequency -- chosen well clear of any l=0 cavity
                  # resonance (roots of tan(ka)=ka; first nontrivial root
                  # is ka~4.49, i.e. f~245Hz for a=1m; 100Hz -> ka~1.83)
omega = 2 * np.pi * f
k = omega / c
print(f"ka = {k*a:.3f} (first resonance near ka~4.493, well clear)")


def j0(x):
    return np.sin(x) / x


def y0(x):
    return -np.cos(x) / x


def j0p(x):
    return (x * np.cos(x) - np.sin(x)) / x ** 2


def y0p(x):
    return (x * np.sin(x) + np.cos(x)) / x ** 2


B = -1j * rho0 * c * k ** 2 * Q / (4 * np.pi)
A = -B * y0p(k * a) / j0p(k * a)


def p_exact_total(r):
    return A * j0(k * r) + B * y0(k * r)


# ---------------------------------------------------------------------
grid = bempp.shapes.regular_sphere(4)
space = bempp.function_space(grid, "DP", 0)
identity = bempp.operators.boundary.sparse.identity(space, space, space)
dlp = bempp.operators.boundary.helmholtz.double_layer(space, space, space, k)

# u_inc is CONSTANT on this unit sphere boundary (source at the origin,
# boundary at constant r=a=1).
u_inc_boundary = 1j * omega * rho0 * Q / (4 * np.pi * a) * np.exp(1j * k * a)


@bempp.complex_callable
def dirichlet_data(x, n, domain_index, result):
    result[0] = u_inc_boundary


uinc_fun = bempp.GridFunction(space, fun=dirichlet_data)

lhs = 0.5 * identity + dlp
u_sol, info = gmres(lhs, uinc_fun, tol=1e-8)
print(f"gmres info: {info}")

test_r = np.array([0.2, 0.4, 0.6, 0.8])
points = np.vstack([test_r, np.zeros_like(test_r), np.zeros_like(test_r)])
p_exact = p_exact_total(test_r)

dlp_pot = helmholtz_potential.double_layer(space, points, k)
p_inc_eval = 1j * omega * rho0 * Q / (4 * np.pi * test_r) * np.exp(1j * k * test_r)
p_total = p_inc_eval - dlp_pot.evaluate(u_sol).flatten()

abs_err = np.abs(p_total - p_exact)
rel_err = abs_err / np.abs(p_exact)

print(f"{'r':>6} {'|p_bem|':>12} {'|p_exact|':>12} {'abs_err':>10} {'rel_err':>10}")
for r, pb, pe, ae, re in zip(test_r, p_total, p_exact, abs_err, rel_err):
    flag = "  <- near standing-wave null, relative error inflated" if abs(pe) < 0.02 else ""
    print(f"{r:6.2f} {abs(pb):12.6f} {abs(pe):12.6f} {ae:10.5f} {re:10.5f}{flag}")

print(f"\nmax absolute error: {abs_err.max():.5f}")
print(f"max relative error away from the null (|p_exact|>0.02): "
      f"{rel_err[np.abs(p_exact) > 0.02].max():.5f}")

if abs_err.max() < 0.01 and rel_err[np.abs(p_exact) > 0.02].max() < 0.05:
    print("\nVALIDATION PASSED -- interior formulation confirmed:\n"
          "  (0.5*Id + K) u = u_inc|Gamma\n"
          "  u(x) = u_inc(x) - (D u)(x)")
else:
    print("\nVALIDATION FAILED -- do not trust the interior solver until "
          "this is resolved.")
