"""
Shared constants for the v2 (Clondalkin-specific) Round Tower rebuild.

Geometry source: RIA archive cross-section drawing of Clondalkin round
tower (Ireland's narrowest surviving round tower), transcribed into
rebuild_spec.md. This module is the single source of
truth for v2 geometry/sweep parameters -- v1's scripts (build_tower_mesh.py,
solve_interior_field.py, solve_tower_radiation.py) each hardcoded their own
copies of shared constants (a real drift risk v1 lived with); v2 avoids
that by having every v2_*.py / *_v2.py script import from here.

SCOPE THIS ROUND: Variant A (ceiled chamber, rigid flat ceiling at the
corbel line) only, at the single nominal wall thickness (0.745m) only.
Variant B (open roof, chamber continuous with the conical void) and the
wall-thickness sensitivity cases (0.64m / 0.813m) are deferred to a
follow-up round -- see rebuild_spec.md SS2/SS5. The
WALL_THICKNESS/VARIANT constants below are still exposed as the intended
single knob for that follow-up work, so build_tower_mesh_v2.py can accept
them as CLI overrides without a rewrite.

v1's round-1 model (superseded by this rebuild -- see
../phase1_generic_tower/PROJECT_SUMMARY_FOR_REVIEW.md) is NOT modified by anything here.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import math

# ---------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------
C_SOUND = 343.0     # m/s
RHO0 = 1.225        # kg/m^3

# ---------------------------------------------------------------------
# Geometry -- RIA drawing, Clondalkin bell storey (SS1 of the instruction)
# ---------------------------------------------------------------------
VARIANT = "A"              # "A" (ceiled, primary) or "B" (open) -- B not
                            # implemented this round, flag kept for the
                            # follow-up round
WALL_THICKNESS = 0.745     # m, derived from two independent measurements
                            # (published base circumference -> 4.043m
                            # external diameter, Devenish-comparable taper
                            # estimate -> 3.472m external diameter at bell
                            # storey). Alternate estimates 0.64m / 0.813m
                            # are explicit sensitivity cases, deferred.

R_INNER = 0.9905           # m, = 1.981m internal diameter / 2 (RIA drawing)
R_OUTER = R_INNER + WALL_THICKNESS   # m, 1.7355 -- bell-storey external
                                      # radius, NOT the tower's wider base
                                      # radius (2.0215m from the 4.043m
                                      # base diameter) -- we only model the
                                      # bell storey, so the base taper is
                                      # irrelevant here.

CHAMBER_HEIGHT = 2.896     # m, floor to ceiling (Variant A) -- RIA drawing
CONE_RISE = 1.905          # m, cone height above the ceiling/springing --
                            # exterior shape only in Variant A (the ceiling
                            # seals the interior domain off from the cone
                            # void; the cone remains part of the real stone
                            # exterior the outdoor air sees regardless)

NUM_WINDOWS = 4
WINDOW_WIDTH = 0.610       # m
WINDOW_HEIGHT = 1.092      # m
WINDOW_SILL_Z = 1.219      # m above floor -- ESTIMATED, not read directly
                            # off the drawing (SS5 item 2 of the instruction);
                            # falls inside the window band unlike v1's
                            # source position, which sat above it
WINDOW_HEAD_Z = WINDOW_SILL_Z + WINDOW_HEIGHT   # 2.311 m

GMRES_MAXITER = 400   # explicit cap, added after investigating why the
                        # first exterior production sweep attempt ran 10+
                        # hours (vs an original ~3-3.5h estimate): a
                        # diagnostic single-solve run (150Hz, with bempp's
                        # logging enabled) showed a HEALTHY solve --
                        # converged in 69 GMRES iterations, ~190s total
                        # (dominated by dense-matrix assembly, not the
                        # iterative solve itself, which took ~9s). But the
                        # SAME production interior sweep already showed
                        # gmres info=31800 at 440/460/480Hz -- i.e. GMRES
                        # failing to converge after roughly 10x the matrix
                        # dimension in iterations, exactly the frequencies
                        # where the mesh's own elem/wavelength margin was
                        # already flagged as marginal (see MESH_SIZE's
                        # comment). Without a cap, a single frequency like
                        # that can run for hours before scipy's default
                        # maxiter (~10*N) gives up -- multiplied across a
                        # 45-point sweep, a few such frequencies plausibly
                        # explain the 10+ hour runaway. 400 is generous
                        # (~6x the healthy 69-iteration case) while bounding
                        # any single bad frequency's cost to roughly the
                        # same order as a normal one; a frequency that
                        # still doesn't converge within that gets flagged
                        # (info!=0, same as before) and the sweep moves on,
                        # rather than silently grinding for hours.

FLOOR_AREA = math.pi * R_INNER ** 2   # 3.084 m^2 -- matches the
                                        # instruction's independently-stated
                                        # 3.08 m^2 (cross-check, see
                                        # PROJECT_SUMMARY_ADDENDUM / plan)

# Window subtense at each face (SS1 of the instruction) -- reported, not
# used directly in the mesh/solve, but kept here so every script that
# wants to print/log it uses the same numbers.
WINDOW_SUBTENSE_EXTERNAL_DEG = math.degrees(WINDOW_WIDTH / R_OUTER)
WINDOW_SUBTENSE_INTERNAL_DEG = math.degrees(WINDOW_WIDTH / R_INNER)
WINDOW_FRACTION_EXTERNAL = NUM_WINDOWS * WINDOW_SUBTENSE_EXTERNAL_DEG / 360.0
WINDOW_FRACTION_INTERNAL = NUM_WINDOWS * WINDOW_SUBTENSE_INTERNAL_DEG / 360.0

# ---------------------------------------------------------------------
# Source (SS3 of the instruction)
# ---------------------------------------------------------------------
SOURCE_POSITION = (0.0, 0.0, 1.50)   # m, chamber axis, standing adult
                                       # mouth height -- now INSIDE the
                                       # window band (1.219-2.311m), unlike
                                       # v1's source, which sat above the
                                       # window tops. One occupant only
                                       # (floor area 3.08 m^2 -- no
                                       # paired/multiple sources modelled).

TARGET_SPL_AT_1M_DB = 84.0   # dB, unchanged from v1 -- not addressed by
                               # the v2 instruction, so not revised here
REFERENCE_FREQ_HZ = 150.0
LOSS_DELTA = 0.01             # k -> k*(1 - i*delta), same placeholder
                               # damping v1 used (still uncalibrated, no
                               # measured reverberation data for this site)

# ---------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------
MESH_SIZE = 0.15   # m -- FALLBACK VALUE, used per the plan's own explicit
                    # contingency ("fall back to v1's 0.15m if 0.12m proves
                    # impractically slow"). Empirical timing on the real v2
                    # exterior mesh at the originally-targeted 0.12m: a
                    # SINGLE solve at 150Hz took 413.6s (12,092 elements) --
                    # far above v1's own ~50s/solve at a similar element
                    # count (v1's exterior solve also assembles two dense
                    # operators, dlp+slp, same as this one, so that's not
                    # the explanation; most likely this environment is
                    # simply slower per-element than v1's original session).
                    # At that rate the ~45-point production sweep would
                    # take 5+ hours for the exterior mesh alone -- not
                    # practical. Reverting to v1's own 0.15m gives
                    # 343/500/0.15=4.57 elem/wavelength at the new 500Hz
                    # ceiling -- BELOW v1's own "6-10 elements/wavelength"
                    # target and below the 5.7-5.87 margin v1 called
                    # "already marginal" at their 390Hz ceiling. This is a
                    # real, explicitly-flagged accuracy compromise, forced
                    # by compute reality rather than silently accepted --
                    # treat gain/SPL figures above ~420-450Hz (where
                    # elem/wavelength drops below ~5) with extra caution,
                    # and prefer the validation gate's own resonance-match
                    # quality as the actual arbiter of trust, not this
                    # ratio alone.

# ---------------------------------------------------------------------
# Analytic validation targets (SS4 of the instruction) -- for
# v2_validate_interior_modes.py. Kept here so the production frequency
# sweep below can be built to include them (see FREQUENCIES).
# ---------------------------------------------------------------------
TRANSVERSE_MODES_HZ = {
    "(1,0)": 101.5,
    "(2,0)": 168.3,
    "(0,1)": 211.2,
    "(3,0)": 231.5,
    "(4,0)": 293.0,
}
AXIAL_MODES_VARIANT_A_HZ = {
    "axial n=1": 59.2,
    "axial n=2": 118.5,
    "axial n=3": 177.7,
    "axial n=4": 236.9,
    "axial n=5": 296.1,
}
# Variant B's quarter-wave axial series (29.6/88.8/148.1 Hz, using
# H=CHAMBER_HEIGHT=2.896m as the "closed" reference plane -- confirmed by
# direct calculation that this, not H=4.801m to the true apex, reproduces
# the given values) is NOT exercised this round (Variant B deferred), kept
# out of FREQUENCIES accordingly, but documented here for the follow-up
# round:
AXIAL_MODES_VARIANT_B_QUARTERWAVE_HZ = {
    "axial n=1": 29.6,
    "axial n=2": 88.8,
    "axial n=3": 148.1,
}

# ---------------------------------------------------------------------
# Production frequency sweep
# ---------------------------------------------------------------------
_V1_FREQUENCIES = [50, 75, 100, 125, 143, 150, 160, 165, 170, 175, 180, 185,
                    190, 196, 200, 207, 210, 220, 230, 240, 250, 260, 270,
                    280, 290, 300, 310, 320, 330, 340, 349, 360, 375, 390]


def _build_frequency_list():
    """v1's 34-point/50-390Hz list, plus every analytic validation target
    (rounded to the nearest Hz) not already within 3Hz of an existing
    point -- so the validation gate's peaks land on an actually-sampled
    point in the production sweep too -- plus evenly-spaced fill points
    every ~20Hz from 400-500Hz (the instruction's new "at least 500Hz"
    ceiling), sorted and de-duplicated (points within 1Hz merged)."""
    points = list(_V1_FREQUENCIES)

    targets = list(TRANSVERSE_MODES_HZ.values()) + list(AXIAL_MODES_VARIANT_A_HZ.values())
    for t in targets:
        t_round = round(t)
        if not any(abs(t_round - p) < 3 for p in points):
            points.append(t_round)

    points += list(range(400, 501, 20))   # 400,420,440,460,480,500

    points.sort()
    deduped = []
    for p in points:
        if not deduped or abs(p - deduped[-1]) >= 1:
            deduped.append(p)
        # else: within 1Hz of the previous point, drop as a duplicate
    return deduped


FREQUENCIES = _build_frequency_list()

# ---------------------------------------------------------------------
# Site / evaluation geometry (unchanged from v1 unless noted)
# ---------------------------------------------------------------------
TOWER_ELEVATION = 30.0     # m, upper chamber floor above ground -- same
                             # 90-110ft-estimate midpoint v1 used; the RIA
                             # drawing didn't address overall tower height
LISTENER_HEIGHT = 1.6       # m

GROUND_DISTANCES = [15, 25, 30, 50, 100, 200, 500, 1000, 2000]   # m --
    # v1's grid (50/100/200/500/1000/2000) plus the 3 real tower-to-church
    # distances from round 1's Herity 1984 addendum (Clonmacnois 30m,
    # Kells 25m, Glendalough 15m), per instruction SS6 step 3

AZIMUTHS_DEG = list(range(0, 360, 10))   # 36 points, 10deg spacing, from
    # the START this time -- round 1 only reached this resolution after
    # discovering its original 8-point/45deg grid aliased with the
    # tower's 4-fold window symmetry and produced a spurious "sweet spot"
    # that had to be retracted. Do not regress to a coarser default here.

ELEVATION_RADIUS = 200.0
ELEVATION_ANGLES_DEG = list(range(-90, 91, 15))
ELEVATION_AZIMUTHS_DEG = [i * 22.5 for i in range(16)]   # 16 points,
    # 22.5deg spacing -- widened from v1's deliberately-coarse 8-point/
    # 45deg elevation grid (a moderate improvement, not the full 36-point
    # resolution used for the ground grid, per the instruction's "use
    # judgement" on this specific sweep)

TOWER_OFFSET = 17.0
RING_RADII = {"2nd_ring": 45.0, "3rd_ring": 91.5}
N_POINTS_AROUND_RING = 36

# ---------------------------------------------------------------------
# Domain tags (physical group numbering) -- same convention v1 used
# ---------------------------------------------------------------------
WINDOWS_DOMAIN = 1
WALLS_DOMAIN = 2

# ---------------------------------------------------------------------
# File paths
# ---------------------------------------------------------------------
MODEL_DIR = f"{HERE}"
MESH_INTERIOR_A_PATH = f"{MODEL_DIR}/v2_mesh_interior_A.msh"
MESH_EXTERIOR_PATH = f"{MODEL_DIR}/v2_mesh_exterior.msh"
DRIVING_CSV_PATH = f"{MODEL_DIR}/v2_window_driving_conditions.csv"


def provenance_header(script_name, extra=None):
    """Returns a list of '# key: value' strings to prepend to every v2
    output file (CSV/plot), per the instruction's SS7 output convention:
    every file states geometry version, mesh variant, wall thickness,
    window width, and producing script."""
    lines = [
        "# geometry_version: v2 (Clondalkin RIA drawing)",
        f"# variant: {VARIANT}",
        f"# wall_thickness_m: {WALL_THICKNESS}",
        f"# window_width_m: {WINDOW_WIDTH}",
        f"# producing_script: {script_name}",
    ]
    if extra:
        # file paths are recorded by name only, so outputs carry no machine-specific paths
        lines.extend(f"# {k}: {os.path.basename(v) if isinstance(v, str) and os.sep in v else v}"
                     for k, v in extra.items())
    return lines


if __name__ == "__main__":
    print(f"R_INNER={R_INNER}  R_OUTER={R_OUTER}  FLOOR_AREA={FLOOR_AREA:.4f} m^2 "
          f"(instruction states 3.08 m^2)")
    print(f"window subtense: external {WINDOW_SUBTENSE_EXTERNAL_DEG:.1f} deg/window "
          f"({WINDOW_FRACTION_EXTERNAL:.1%} of circumference), "
          f"internal {WINDOW_SUBTENSE_INTERNAL_DEG:.1f} deg/window "
          f"({WINDOW_FRACTION_INTERNAL:.1%})  "
          f"(instruction states 20.2deg/22% external, 35.8deg/40% internal)")
    print(f"FREQUENCIES ({len(FREQUENCIES)} points): {FREQUENCIES}")

# Frequencies whose exterior solve hit GMRES_MAXITER without reaching tolerance in the
# v2 sweep (see v2_headline_readout.md section 1). Exclude them from headline numbers.
NONCONVERGED_HZ = {250, 260, 280, 290, 300, 310, 320, 340, 349, 360, 375, 390,
                   400, 420, 440, 460, 480, 500}
