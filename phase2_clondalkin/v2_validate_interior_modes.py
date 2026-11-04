"""
Validation gate for the v2 interior-A mesh -- per
rebuild_spec.md SS4: "The BEM interior solution should
show features near these [analytic rigid-cylinder modes]... if it
doesn't, stop and debug before proceeding."

This is a STANDALONE gate script, meant to be run and read BEFORE
launching the production frequency sweep (solve_interior_field_v2.py's
__main__ block) -- not invoked automatically inside the pipeline.

Method:
  1. Cross-check the instruction's 10 hardcoded analytic target
     frequencies against independently-computed closed-form values
     (scipy Bessel-derivative zeros for the transverse modes, c/(2H)
     series for the axial modes) -- a self-check that the transcribed
     numbers in the instruction file weren't themselves wrong.
  2. Solve the interior-A BEM problem across a uniform, moderately dense
     frequency grid spanning all 10 targets, evaluate the pressure at a
     FIXED off-axis, off-source probe point (interior representation
     u(x) = u_inc(x) - D[u](x), x strictly interior -- same formulation
     validated in v1's validate_interior_bem.py), and peak-pick |u(f)|.
  3. Match each detected peak to its nearest analytic target within a
     tolerance band, write v2_validation.md (pass/fail table) and a plot.

The probe point is deliberately NOT the source location (to avoid
sampling exactly at a mode's own node) and NOT on the z-axis (transverse
modes with azimuthal order m>0 are node lines on-axis).

REVISION after first run: the first attempt used the production model's
actual on-axis source (P.SOURCE_POSITION) and a probe at 45deg azimuth.
8/10 targets passed cleanly (including a PERFECT 0.0Hz match on the
(4,0) mode and the entire axial series) but (1,0) and (2,0) failed. This
is NOT a geometry or solver defect -- it's a real symmetry consequence:
an on-axis source in a cavity with 4-fold (C4) window symmetry can only
couple to axisymmetric (m=0) and 4-fold-compatible (m=4, m=8, ...) modes;
m=1 and m=2 are structurally under-excited by symmetry, near-independent
of probe choice. The 45deg probe additionally sat on an exact NODE for
the m=2 mode (cos(2*45deg)=0), compounding the miss for (2,0) specifically.
Fixed here by using a deliberately OFF-AXIS source for THIS DIAGNOSTIC
CHECK ONLY (VALIDATION_SOURCE_POSITION below) -- the production model's
actual on-axis source (SIF.solve_case / P.SOURCE_POSITION) is completely
untouched by this change. Verifying "does the BEM code reproduce this
cylinder's true resonant structure" is a property of the geometry alone
and is best tested with a generic excitation that couples to every mode
order; separately confirming what a real on-axis singer's source actually
excites is the production interior sweep's job (and is exactly what
SS4.1-style findings in this project's own history are about), not this
gate's.
"""
import numpy as np
import scipy.special
from scipy.signal import find_peaks
import bempp_cl.api as bempp
from bempp_cl.api.linalg import gmres
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import v2_params as P
import solve_interior_field_v2 as SIF   # reuses grid/space/identity/Q setup
                                          # (module-level code only --
                                          # the production __main__ sweep
                                          # is guarded and does NOT run
                                          # on import)

VALIDATION_SOURCE_POSITION = np.array([0.5 * np.cos(np.radians(20)), 0.5 * np.sin(np.radians(20)), 1.3])
# off-axis (r=0.5m, well inside R_INNER), 20deg azimuth -- deliberately
# NOT a multiple of 90deg (window bearing), 45deg (a node for m=2), or
# any other "special" angle -- so this source couples to every azimuthal
# mode order, unlike the production model's on-axis source (see docstring
# revision note above).

PROBE_POINT = np.array([0.6 * np.cos(np.radians(200)), 0.6 * np.sin(np.radians(200)), 0.8])
# r=0.6m, 200deg azimuth (roughly opposite the source, also off any
# special angle), z=0.8m (different from both the validation source's
# z=1.3m and the production source's z=1.5m).

FREQ_MIN_HZ = 45
FREQ_MAX_HZ = 310
FREQ_STEP_HZ = 4.0   # covers all 10 targets' neighbourhoods. Widened from
                      # an initial 3Hz plan after empirical timing on the
                      # real v2 meshes came in well above estimate (see
                      # v2_params.MESH_SIZE's comment) -- 4Hz still safely
                      # resolves peaks broadened by LOSS_DELTA=0.01 damping
                      # (expected half-width on the order of a few Hz at
                      # these frequencies) while cutting point count ~25%.

TOLERANCE_FRAC = 0.05
TOLERANCE_MIN_HZ = 8.0

OUT_MD = f"{P.MODEL_DIR}/v2_validation.md"
OUT_PLOT = f"{P.MODEL_DIR}/v2_validation_modes.png"

# ---------------------------------------------------------------------
# Step 1: independent analytic cross-check
# ---------------------------------------------------------------------
_TRANSVERSE_M = {"(1,0)": 1, "(2,0)": 2, "(0,1)": 0, "(3,0)": 3, "(4,0)": 4}


def crosscheck_transverse():
    print("--- transverse mode cross-check (scipy Bessel-derivative zeros) ---")
    ok = True
    for label, f_given in P.TRANSVERSE_MODES_HZ.items():
        m = _TRANSVERSE_M[label]
        zero = scipy.special.jnp_zeros(m, 1)[0]
        f_calc = P.C_SOUND * zero / (2 * np.pi * P.R_INNER)
        delta = abs(f_calc - f_given)
        status = "OK" if delta < 0.5 else "MISMATCH"
        if delta >= 0.5:
            ok = False
        print(f"  {label}: given={f_given:.1f}Hz  calc={f_calc:.2f}Hz  "
              f"delta={delta:.2f}Hz  [{status}]")
    return ok


def crosscheck_axial():
    print("--- axial mode cross-check (closed-closed cylinder, f_n=n*c/(2H)) ---")
    ok = True
    for label, f_given in P.AXIAL_MODES_VARIANT_A_HZ.items():
        n = int(label.split("=")[1])
        f_calc = n * P.C_SOUND / (2 * P.CHAMBER_HEIGHT)
        delta = abs(f_calc - f_given)
        status = "OK" if delta < 0.5 else "MISMATCH"
        if delta >= 0.5:
            ok = False
        print(f"  {label}: given={f_given:.1f}Hz  calc={f_calc:.2f}Hz  "
              f"delta={delta:.2f}Hz  [{status}]")
    return ok


# ---------------------------------------------------------------------
# Step 2: dense probe sweep
# ---------------------------------------------------------------------
def solve_and_probe(f, probe_pt):
    omega = 2 * np.pi * f
    k_real = omega / P.C_SOUND
    k = k_real * (1 - 1j * P.LOSS_DELTA)

    dlp = bempp.operators.boundary.helmholtz.double_layer(SIF.space, SIF.space, SIF.space, k)
    lhs_op = 0.5 * SIF.identity + dlp
    src = VALIDATION_SOURCE_POSITION   # off-axis, deliberately NOT
                                         # P.SOURCE_POSITION -- see module
                                         # docstring's revision note

    @bempp.complex_callable
    def dirichlet_data(x, n, domain_index, result):
        dx = x - src
        r = max(np.linalg.norm(dx), 1e-6)
        result[0] = 1j * omega * P.RHO0 * SIF.Q / (4 * np.pi * r) * np.exp(1j * k * r)

    uinc_fun = bempp.GridFunction(SIF.space, fun=dirichlet_data)
    u_sol, info = gmres(lhs_op, uinc_fun, tol=1e-5)
    if info != 0:
        print(f"  [WARNING] f={f}Hz probe solve: gmres info={info}")

    pt = np.array(probe_pt).reshape(3, 1)
    dlp_pot = helmholtz_potential.double_layer(SIF.space, pt, k)
    Du = dlp_pot.evaluate(u_sol).flatten()[0]

    dx = probe_pt - src
    r = max(np.linalg.norm(dx), 1e-6)
    u_inc = 1j * omega * P.RHO0 * SIF.Q / (4 * np.pi * r) * np.exp(1j * k * r)
    return u_inc - Du


def run_probe_sweep():
    freqs = np.arange(FREQ_MIN_HZ, FREQ_MAX_HZ + FREQ_STEP_HZ / 2, FREQ_STEP_HZ)
    mags = []
    for f in freqs:
        u = solve_and_probe(f, PROBE_POINT)
        mags.append(abs(u))
        print(f"  probe f={f:.1f}Hz  |u|={abs(u):.6f}")
    return freqs, np.array(mags)


# ---------------------------------------------------------------------
# Step 3: peak matching
# ---------------------------------------------------------------------
def match_targets(freqs, mags, all_targets):
    peak_idx, _ = find_peaks(mags)
    peak_freqs = freqs[peak_idx]

    results = []
    for label, f_target in all_targets:
        if len(peak_freqs) == 0:
            nearest = None
            delta = None
        else:
            i = np.argmin(np.abs(peak_freqs - f_target))
            nearest = peak_freqs[i]
            delta = abs(nearest - f_target)
        tol = max(TOLERANCE_MIN_HZ, TOLERANCE_FRAC * f_target)
        passed = (delta is not None) and (delta <= tol)
        results.append((label, f_target, nearest, delta, tol, passed))
    return results


def write_report(results, analytic_ok, freqs, mags):
    with open(OUT_MD, "w") as fh:
        fh.write("# v2 validation gate -- interior-A analytic mode comparison\n\n")
        fh.write(f"Probe point: {PROBE_POINT.tolist()} (off-axis, off-source)\n\n")
        fh.write(f"Analytic-target self-check (transcribed vs independently computed): "
                 f"{'PASS' if analytic_ok else 'FAIL -- see console log'}\n\n")
        fh.write("| Target | Target Hz | Nearest detected peak Hz | Delta Hz | Tolerance Hz | Result |\n")
        fh.write("|---|---|---|---|---|---|\n")
        all_pass = True
        for label, f_target, nearest, delta, tol, passed in results:
            all_pass = all_pass and passed
            nearest_s = f"{nearest:.1f}" if nearest is not None else "none detected"
            delta_s = f"{delta:.1f}" if delta is not None else "-"
            fh.write(f"| {label} | {f_target:.1f} | {nearest_s} | {delta_s} | {tol:.1f} | "
                     f"{'PASS' if passed else 'FAIL'} |\n")
        overall = "PASS" if all_pass else ("FAIL -- stop and debug before trusting "
                                            "the production sweep, per the instruction")
        fh.write(f"\n**Overall: {overall}**\n")
        fh.write(f"\nProbe sweep: {freqs[0]:.0f}-{freqs[-1]:.0f}Hz, step {FREQ_STEP_HZ}Hz, "
                 f"{len(freqs)} points, mesh v2_mesh_interior_A.msh (r_inner={P.R_INNER}m).\n")
    return all_pass


def plot_report(freqs, mags, all_targets, out_path):
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(freqs, mags, "-", color="#2a4d69", lw=1.3, label="|u(probe)|")
    for label, f_target in all_targets:
        ax.axvline(f_target, color="#c0392b", ls="--", lw=0.8, alpha=0.7)
        ax.annotate(label, (f_target, ax.get_ylim()[1] * 0.95), rotation=90,
                    fontsize=7, color="#c0392b", ha="right", va="top")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("|pressure| at probe point (arb. units)")
    ax.set_title("v2 validation gate: interior-A probe response vs analytic rigid-cylinder modes")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    ok1 = crosscheck_transverse()
    ok2 = crosscheck_axial()
    analytic_ok = ok1 and ok2

    all_targets = (list(P.TRANSVERSE_MODES_HZ.items())
                   + list(P.AXIAL_MODES_VARIANT_A_HZ.items()))
    all_targets.sort(key=lambda t: t[1])

    print(f"\n--- probe sweep: {FREQ_MIN_HZ}-{FREQ_MAX_HZ}Hz, step {FREQ_STEP_HZ}Hz ---")
    freqs, mags = run_probe_sweep()

    results = match_targets(freqs, mags, all_targets)
    all_pass = write_report(results, analytic_ok, freqs, mags)
    plot_report(freqs, mags, all_targets, OUT_PLOT)

    print(f"\n=== VALIDATION GATE: {'PASS' if all_pass else 'FAIL'} ===")
    print(f"see {OUT_MD}")
    if not all_pass:
        raise SystemExit(1)
