# Round Tower Acoustic Modelling — Summary for Independent Review

## 0. MAJOR CORRECTION — read this before anything else in this document

A follow-up review pass (`review_round_notes.md`) flagged that the
generic distance/azimuth evaluation grid used 8 azimuth samples at 45°
spacing — an exact submultiple of the tower's 90° (4-fold) window symmetry,
which risks aliasing with the real angular lobe structure rather than
randomly missing it. This was checked and confirmed real: **widening to 36
azimuth samples (10° spacing) substantially flattened the frequency-response
curve that most of §4 below was built on.**

What this changes: the earlier claimed sharp peak near 185–196Hz (gain
approaching 0dB) and sharp trough near 250Hz (gain ≈ −12 to −15dB) — the
entire basis for the Martyrs/French note-overlay narrative ("tonic sits on
the peak, thirds sit in the trough") — **shrinks to a much flatter curve**
ranging roughly −5.2dB (best, now at 125Hz) to −9.8dB (worst, at 50Hz, the
low end, not 250Hz). At 196Hz specifically, properly-resolved gain is
−7.7dB — one of the *weaker* points in the mid-range, not the best. The
Martyrs/French tonic-vs-third contrasts shrink from ~10dB to ~1–2dB.

**What survives:** a small, consistent tonic-favouring tendency (the tonic
of both tunes checked still comes out slightly better than their thirds,
just not dramatically so), the near-perfect 4-window balance from a centred
source (§4.1, unaffected — that's an interior-field finding, not dependent
on this exterior evaluation grid), and the qualitative existence of sharp
directional interference structure (confirmed real, see §4.7) — just not
the specific magnitude of peak/trough claimed earlier.

**How to read the rest of this document:** §4.2–§4.4 and §4.6 below were
written before this correction and contain numbers that are now superseded.
They are left in place (not silently edited) so the revision is visible,
each marked SUPERSEDED where it occurs. §4.7 and §4.8 (new) contain the
corrected numbers and the additional checks (reciprocity, CFIE cross-check)
that came out of the same review pass. If you only read one section for
current, trustworthy numbers, read §4.7. §4.9 (new) is a follow-on: once
§4.7 removed the sharp resonance story, the project's focus shifted from
long-range projection toward a narrower, more practical question — does
the tower's height help a voice be heard over a crowd gathered nearby —
covered there with its own, separately-flagged uncertainty.

## 1. Hypothesis under test

Irish Round Towers (9th–12th c. CE) are conventionally explained as bell towers,
but this is disputed (no bells large enough have been found in association with
them; large cast bells arrive in Ireland in the 12th century, coincident with
the towers going out of use). The alternative hypothesis under test here: the
upper chamber of a Round Tower — a cylindrical room with a corbelled conical
roof and four apertures at the cardinal compass points — functioned as a
passive acoustic device, converting a single standing human voice (a lector or
precentor chanting) into an efficient, roughly even, four-directional outdoor
broadcast, audible to a congregation gathered outside the tower's base.

**Revision note**: an earlier version of this model used a shallow-DOME roof
(a guess, no source). That guess was replaced with a straight-sided CONE
after the user obtained an actual (if imperfect) 3D scan measurement of
Clondalkin round tower's intact interior roof. This is a substantial
geometry change — cone vs. dome are meaningfully different acoustic
reflectors — and re-running the full pipeline on the corrected geometry
was treated as a deliberate robustness check on the earlier findings (see
§4.6). The dome-based results are superseded; where relevant this document
notes what changed and what didn't.

This is a computational acoustics investigation (Boundary Element Method) of
that specific mechanical/acoustic claim. It does not address the historical or
liturgical plausibility of the hypothesis independently — only whether the
geometry, as far as it is known, behaves the way the hypothesis requires.

## 2. What is measured/real vs. what is assumed — read this first

This distinction matters more than any single number below.

**Measured / survey-derived (not guessed):**
- Roof shape: a straight-sided CONE ("regular cone of smooth corbels"),
  from an actual (imperfect) 3D scan of Clondalkin round tower's intact
  interior roof, replacing an earlier shallow-dome guess
- Roof/window geometry, Clondalkin scan: cone springing sits just above the
  window tops (~zero blank wall between window top and roof start); apex
  8ft6in (2.591m) above the window tops; window height ~4ft (1.219m,
  stated directly) though NOTE a real internal inconsistency — the scan's
  two apex-referenced distances (8ft6in to window top, 12ft to window
  bottom) subtract to 3ft6in (1.067m), ~0.15m shorter than the direct 4ft
  estimate; the direct estimate was used, not silently resolved; window
  width ~2ft (0.610m, "about")
- Windows centred on the four cardinal compass points, bases at floor level
- Site geometry: three concentric monastic enclosure rings at radii 20m /
  45m / 91.5m from a common centre; tower positioned at the inner edge of the
  2nd ring, 17m west of that centre, tower door (and by inference window 0)
  facing east toward the church. **Caveat added after real comparative data
  became available (see §4.9's Herity 1984 addendum): this "due west" bearing
  is this project's own chosen convention for its one modeled site, not a
  general rule.** Real sites vary in *direction*, not just distance — the
  three round-tower sites with quantified positions in Herity (1984) are all
  offset NW or SW of the church (15–30m), never exactly due west, though all
  keep the tower door targeting the church across a shared open space. This
  project's exact-west convention was a simplification adopted before that
  data was in hand; the underlying modeled geometry wasn't re-run to match,
  but the assumption should now be read as "one reasonable case," not a
  documented fact about round towers generally.
- Tower's own footprint: ~4.5–5m diameter, matching the model's 2.2m radius
  closely (this radius figure is NOT from the Clondalkin scan specifically
  — see below)
- Gaelic psalm precentor vocal range and Sean-nós solo male vocal range
  (given by the user, presented as informed estimates, not the user's own
  primary-source citations — treat as secondary data of unknown provenance)

**Explicitly guessed / assumed, clearly flagged in code comments throughout:**
- Interior chamber radius: 2.2m — a compromise between the tower's inner
  wall radius (1.9m) and outer radius (2.5m) from an earlier, different,
  general-survey-average source (not Clondalkin-specific), used because the
  acoustic model is a zero-thickness shell (wall thickness itself is not
  modelled — argued to be immaterial for exterior radiation, see §3). Mixed
  provenance, worth flagging: the roof and window dimensions now come from
  one specific tower's scan (Clondalkin) while the radius still comes from
  a general average across towers — not fully self-consistent as "one
  real tower's measurements."
- Source acoustic power: calibrated to a target 84dB SPL at 1m (revised down
  from an initial 95dB guess, which sat in "exceptional peaks" territory
  per published SPL-at-1m tables; 84dB sits between "amateur singing"
  (70-80dB) and "trained classical male singer" (80-90dB), consistent with
  muezzin call-to-prayer measurements ~80-90dB at 1m — still a guess, no
  measurement exists for sean-nós/precentor delivery specifically, but
  better-grounded), held as a constant volume-velocity monopole across the
  frequency sweep. Because the whole BEM pipeline is linear in source
  amplitude, this revision shifted every absolute SPL figure in this
  document down by exactly 11dB and changed NONE of the relative findings
  (gain vs. free-field, window balance, resonance-band locations) — those
  are all ratios/differences the source-level constant cancels out of.
- Interior acoustic damping: a small uniform loss factor folded into the
  wavenumber (k → k·(1−iδ), δ≈0.01) to keep resonant peaks finite; not
  calibrated against any measured reverberation data (none exists for this
  site)
- Tower elevation above ground: 30m, the midpoint of a 90–110ft estimate
  range given by the user (no independent survey source)

**Known-uncertain by construction:**
- Windows are treated as vibrating/pressure-radiating pistons standing in
  for the "escaping chant," calibrated from a computed interior sound field
  — see §3 for exactly how, and why this is a two-step (not fully
  self-consistent) approximation

## 3. Computational method

### 3.1 Geometry / mesh
Built parametrically in Gmsh (OpenCASCADE kernel) as a closed, zero-thickness
shell: a cylinder (radius 2.2m, height 1.219m) capped with a straight-sided
CONE (apex 2.591m above the springing, tapering to a true point — matching
"smooth corbels up to the apex"), with the four window regions **not cut as
real holes** — they are imprinted as separate physical-surface patches via a
boolean fragment operation. Two physical groups result: "windows" and
"walls" (dimension-2 surfaces in both cases; an earlier hand-built mesh had a
bug where "windows" ended up as dimension-1 curves — a Gmsh GUI classic
error, fixed here by scripting the mesh generation instead of manual
point-and-click). The classification logic that separates "windows" from
"walls/roof" surfaces after the fragment operation was re-derived for the
cone (a linear radius-vs-height taper check) — see `build_tower_mesh.py`,
`classify_skin()`.

An earlier version of this geometry used a shallow dome (radius-scaled
hemisphere, 0.9m rise) instead of the cone — a guess, since no source gave
the roof shape. That version's results are superseded; §4.6 documents the
comparison.

Two mesh resolutions were used at different points in this project:
- 0.2m elements (~5.7 elements/wavelength at 300Hz) — inherited from the
  original (incomplete) methodology document this project started from,
  which itself targeted "6–10 elements per acoustic wavelength" up to
  ~300Hz. This was already marginal at its own ceiling.
- 0.15m elements (~5.7 elements/wavelength at 400Hz) — re-cut later to
  extend the valid frequency range to ~390Hz, to cover higher notes in the
  vocal ranges being compared against. A finer 0.1m mesh was tried first but
  found computationally impractical (see §3.4). The cone geometry's smaller
  total surface area (shorter cylindrical wall, though a taller roof)
  produces a smaller mesh at the same 0.15m target size than the dome
  version did (~4,700 vs. ~6,200 elements).

### 3.2 Exterior radiation — validated formulation
Exterior sound radiation from the tower (Helmholtz equation, unbounded
exterior domain, radiation condition at infinity) is solved via the direct
boundary integral equation:

    (0.5·Id − K) p = −V g

where p is the unknown surface pressure, g = iωρ₀v_n is prescribed Neumann
(velocity) data, K is the double-layer operator, V is the single-layer
operator (both `bempp_cl.api.operators.boundary.helmholtz`). Exterior field
at a point x is then p(x) = D[p](x) − S[g](x) (potential operators).

**This exact formulation was validated before use**, not assumed correct:
checked against the closed-form analytical solution for a pulsating sphere
(uniform prescribed radial velocity on a rigid sphere of known radius),
achieving <0.6% error in both magnitude and phase across 5 test radii and a
representative frequency. See `validate_bem_formulation.py`.

Known limitation of this formulation: being a direct (non-combined/non-CFIE)
method, it can in principle fail at isolated "irregular frequencies"
coinciding with interior Dirichlet eigenvalues of the enclosed volume. Every
solve checks GMRES convergence (`info==0`) and would flag (not silently
accept) any frequency that failed to converge cleanly. No such failures
occurred across the ~34-point frequency sweep used in the final results.

### 3.3 Interior sound field — validated formulation
Rather than assume a boundary condition at the windows (the original
approach), the interior sound field produced by an actual point-source
monopole (standing in for a chanting voice) inside the rigid, closed
chamber is solved directly, and the resulting pressure at each window is
read off and used to drive the exterior problem. This decomposes the field
as p = p_inc + p_scat (p_inc = free-space monopole field of the source
itself), reducing (using the rigid-wall condition ∂p_total/∂n=0 everywhere)
to:

    (0.5·Id + K) u = u_inc|Γ        [u = TOTAL boundary pressure trace]
    u(x) = u_inc(x) − (D u)(x)      [x strictly interior]

**This was also validated before use, against a different closed-form
target**: a monopole source at the exact centre of a rigid sphere excites
only the spherically-symmetric mode, giving a closed-form spherical-Bessel
solution. Validated to <0.3% error away from a standing-wave null (where
relative error is naturally inflated by a small denominator; absolute error
stayed ~2×10⁻⁴ throughout). See `validate_interior_bem.py`. Note: an
earlier attempt at this validation, going through a "scattered field"
formulation with the single-layer operator, produced large (~40%) errors
that could not be resolved by sign-flipping alone; it was abandoned in
favour of the simpler direct-trace formulation above, which validated
cleanly on the first physically-correct attempt.

**This is a weakly-coupled, two-step approximation**, not a single
self-consistent transmission solve: the interior field is computed treating
the windows as closed/rigid (same as the rest of the shell), ignoring the
~5% of shell area they actually occupy and any energy that would leak out
through them and modify the interior field. The resulting pressure is then
converted to a velocity boundary condition (v = p/(ρ₀c), same plane-wave
impedance relation used throughout) and used to drive the exterior problem
as a piston. This is a standard simplification for this class of problem,
reasonable given the window area fraction, but not exact.

### 3.4 Compute cost / practical constraints (relevant to judging thoroughness)
No OpenCL/GPU or FMM (fast multipole) acceleration is available in this
environment (checked; `bempp_cl.api.check_for_fmm()` returns False, and the
`exafmm` package is not available via pip). All BEM assembly is dense,
CPU-only, via Numba JIT. This was the binding constraint on mesh resolution:
a 0.1m mesh (~13,500 elements) made a single-frequency solve take ~300s,
confirmed (by direct timing across 4 different frequencies and 2 different
GMRES tolerances) to be dominated by dense matrix assembly cost, not solver
iteration count or frequency-dependent conditioning. A 0.15m mesh (~6,200
elements) reduced this to ~50s/solve, making a ~34-frequency sweep feasible
(~30–55 minutes per full run). This directly limited how far the frequency
range could be extended and how finely it could be sampled.

## 4. Key results

### 4.1 The room's own symmetry, not assumption, produces even 4-window output
The original (superseded) version of this model simply assumed equal,
in-phase driving at all 4 windows. Once the interior field was actually
solved for a source at the room's centre (~1.75m height, standing, natural
posture — the windows top out at 1.45m, below standing mouth height, so
there is no natural "aim at a window" posture available to a standing
singer), the result **falls out of the physics**: because a source on the
central axis of a near-axisymmetric cavity can only excite rotationally-
symmetric interior modes, the 4 windows come out driven within **0.09dB of
each other on average across 50–390Hz** (worst case 0.21dB at any single
frequency tested).

A control case (source moved near window 0, seated height ~0.7m, reproducing
"sitting and singing straight out one window" — a configuration the user
reports real bell tests showed no broader effect from) was run at the
original coarser mesh/frequency range and showed the opposite: **11.2dB
average imbalance** across windows, with the near window ~4–5× the pressure
of the far windows. (This control case was not re-run at the extended
0.15m-mesh/390Hz range — see §5, open items.)

### 4.2 Exterior radiation is generally somewhat WORSE than an idealised free-field source
Comparing the tower's radiated field to a free-field monopole of the same
total volume velocity (a "fair" baseline — same source strength, no tower
structure), the tower under-performs across most of the tested range. This
does not by itself support a "the tower amplifies the voice" reading — if
anything the solid tower body attenuates relative to an idealised
unobstructed source of equal power, most likely via diffraction/shadowing.
The interesting structure is in *where* it under-performs least (see 4.3).

### 4.3 [SUPERSEDED BY §4.7 — 8-azimuth sampling artifact, see §0] A broad, validated resonance band, found from geometry ALONE, before any vocal data was introduced
[CONE-ROOF GEOMETRY, current model.] The frequency response (gain vs.
free-field baseline, computed purely from the tower's geometry — no
assumption about what pitch anyone actually sings) shows two elevated
bands and one sharp trough:
- **~100–220Hz**: broad peak, best gain ≈ −1.4dB at 185Hz — the single best
  point in the whole sweep
- **~220–260Hz**: sharp trough, worst gain in the model (≈ −12.2dB at 250Hz)
- **~290–330Hz**: second, smaller peak, gain again approaches −3dB

This was confirmed to be a genuine broad feature, not a numerical artifact,
by fine-sampling ~20 closely-spaced points through 100–390Hz — all
converged cleanly (no GMRES warnings), and the elevated region spans many
adjacent sample points smoothly rather than showing an isolated spike.
**This frequency-response shape was computed and characterised before any
vocal pitch data was introduced into the analysis pipeline** — the
comparison in §4.4 was a check against an independently-derived prediction,
not a fit.

**This same computation was done twice, on two materially different roof
geometries (shallow dome, then cone) — see §4.6.** The peak/trough
locations shifted modestly (best point moved from 196Hz/dome to 185Hz/cone;
trough from 240Hz/dome to 250Hz/cone) but the overall two-peak-one-trough
shape, and critically its overlap with the vocal tonic ranges in §4.4,
survived the change. Numbers below are the current (cone) values.

### 4.4 [SUPERSEDED BY §4.7 — see §0] Comparison against vocal tradition pitch data (post-hoc check)
Two sets of vocal range data were provided by the user (their own
characterisation, not independently sourced by this analysis) and overlaid
on the pre-existing frequency-response curve:

**(a) Gaelic psalm-singing (congregational, precentor + crowd, mixed
range)** — subsequently identified by the user as the wrong scope for this
model (the congregation's return is not sound the tower radiates; see
§4.5). Kept here for completeness: two traditional tune profiles (Martyrs,
D Dorian; French, major key) were checked. Both tunes' tonic notes landed
in the model's best-gain zone. Martyrs' "emotional core" note (minor 3rd,
~25% of performance time, heavily held/ornamented) landed in the worst
trough (~230–250Hz) — a genuine tension in the result, not glossed over.
French's structure aligned more favourably: both its dominant pillars
(tonic 35% + fifth "camped out" on, 12.5%, combined ~47.5% of performance
time) landed in good zones; only its major 3rd (20%) fell in the same
trough. Because a minor 3rd (Martyrs) and major 3rd (French) are only one
semitone apart and both tunes are pitched near the same tonic, this "third
lands in the trough" result is understood to be **one underlying finding
(the trough's location relative to a G3/Ab3-ish tonic), not two independent
confirmations** — flagged explicitly to avoid over-counting.

**(b) Sean-nós solo male vocal range (the actually relevant comparison)** —
per the user's clarification, the tower model concerns only a single male
lector/precentor's voice; the congregation's full-range mixed-gender return
happens in the 2nd ring, not from the tower, and is acoustically irrelevant
to what the tower radiates (relevant only as background noise the lector
must be heard over — a separate, unmodelled question, see §4.5). Given
range: tonic C3–G3 (130–196Hz), "sweet spot" D3–E3 (146–165Hz), typical
performance hovering 110–250Hz. Checked against the model:

| Note | Freq | SPL@100m | Gain vs free-field |
|---|---|---|---|
| C3 (low tonic bound) | ~130Hz | ~35dB | ~−1.9dB |
| D3 (sweet-spot start) | 146Hz | ~37dB | ~−2.4dB |
| E3 (sweet-spot end) | 165Hz | 43.5dB | −2.2dB |
| G3 (high tonic bound) | 196Hz | 42.3dB | −1.7dB |
| (model's actual best point) | 185Hz | ~40dB | **−1.4dB (best in the model)** |

(SPL figures shown at 84dB@1m source calibration, CONE geometry; gain vs.
free-field is unchanged by the source-level assumption, see §2. Note the
model's single best-projecting frequency, 185Hz, now falls WITHIN the
stated tonic range rather than right at its upper edge as in the earlier
dome version — if anything a slightly better alignment.)

The entire stated tonic range performs solidly-to-excellently; the top of
that range coincides with the model's single best-performing frequency.

**Checked at the actual pilgrim positions** (not just a generic distance
grid): re-evaluated at points sampled every 10° around the real 2nd ring
(28–62m from the tower, accounting for the tower's 17m off-centre position)
and 3rd ring (74.5–108.5m), across the sean-nós tonic range:

| Ring | Freq range | SPL across the ring (dB) |
|---|---|---|
| 2nd ring (28–62m) | 125–196Hz | 17.5–51.3 depending on position/pitch |
| 2nd ring (28–62m) | at 196Hz (G3) specifically | 39.8 (mean) |
| 3rd ring (74.5–108.5m) | 125–196Hz | 19.7–48.8 depending on position/pitch |
| 3rd ring (74.5–108.5m) | at 196Hz (G3) specifically | 37.7 (mean) |

(CONE geometry, 84dB@1m source calibration. This table was regenerated
efficiently after the geometry change: the pilgrim-ring evaluation was
merged into the same solve pass as the generic distance grid, rather than
re-solving the whole BEM problem a second time purely to evaluate different
points — see `solve_tower_radiation.py`'s `evaluate_pilgrim_rings()`. The
earlier standalone `evaluate_pilgrim_rings.py` script is now superseded by
this merged version.)

Position around the ring matters substantially at a given pitch (near-field
interference between the 4 coherent window sources, not simple distance
falloff — confirmed by checking raw distance-vs-SPL numbers directly rather
than trusting a smooth-falloff assumption), but the tonic range is
consistently in an audible-to-marginally-audible band (roughly quiet-room
level, weaker toward the low end) throughout both rings.

The
2.5–3.5kHz "singer's formant" that gives sean-nós its perceptually piercing/
nasal quality is far above anything this mesh resolves (390Hz ceiling) —
the model can speak to fundamental-pitch projection efficiency, not to the
timbral mechanism that makes the voice cut through.

### 4.5 Elevation/vertical directivity, and what remains unmodelled
Power-weighted, roughly a third of radiated energy goes into the lower
hemisphere (below the tower's own local horizontal) vs. two-thirds upper —
i.e. the model does **not** show a strong deliberate downward bias; it reads
more as "peaked near-horizontal, symmetric rolloff up and down" than
"aimed at the ground." This tempers a simple "obviously aimed at the
congregation" reading.

Separately (raised by the user): whether the tower's ~30m elevation helps
the lector's voice be heard *over* a crowd of hundreds of pilgrims in the
2nd ring, independent of pitch-matching. **This has since been modelled
and revised twice — see §4.9 for the full account.** The intuition stated
in this paragraph originally (crowd = attenuating obstacle, grazing sound
loses energy, elevation avoids that) turns out to be frequency-dependent,
not a fixed "general acoustic principle" as first assumed: at the low end
of this project's frequency range, human bodies are close to acoustically
transparent (wavelength much larger than body size) and the crowd's
effect is small in either direction; only as frequency climbs toward and
through the upper part of the vocal range do bodies become effective
obstacles and the original intuition starts to hold, and hold strongly.
The BEM setup itself still treats the whole propagation path as clear open air (no
crowd representation in the BEM solve); the correction is applied as a
post-processing step, detailed fully in §4.9.

Also explicitly out of scope so far, per the user: proximity to water
(reflective, low-absorption ground surface, common near many Round Tower
sites) as a possible propagation enhancer — flagged for future work, not
modelled.

### 4.6 [NUMBERS SUPERSEDED BY §0/§4.7 -- both dome and cone comparisons below used the flawed 8-azimuth grid; the dome geometry was never re-run at 36-azimuth resolution, so a like-for-like re-comparison is an open item, not done here] Robustness check: shallow dome vs. cone roof
This is treated as a first-class result, not a footnote. The entire
pipeline (interior solve, exterior solve, both note-overlay comparisons)
was originally built and validated using a shallow-dome roof — a pure
guess, since no source specified the roof shape. When the user later
obtained an actual 3D-scan measurement of a different, real tower
(Clondalkin) showing a straight-sided corbelled CONE instead — a
meaningfully different acoustic reflector, springing at a different height,
rising roughly 3x higher above the window band than the dome guess had —
the full pipeline was re-run on the corrected geometry as a deliberate
check on whether the earlier findings were an artifact of the wrong roof
assumption.

**What changed:** low-frequency performance (50–75Hz) got meaningfully
worse (roughly 3–4dB) with the cone; the trough shifted from 240Hz to
250Hz and softened by about 3dB (−15.1dB dome vs. −12.2dB cone); the
best-gain frequency shifted down from 196Hz to 185Hz.

**What did not change:** the overall two-peak/one-trough shape of the
frequency response; the ~100–220Hz broad projection-efficiency band; the
secondary ~290–330Hz peak; the near-perfect 4-window balance from a
centred source (§4.1, confirmed independent of roof shape — expected,
since that result follows from axisymmetry, not roof geometry
specifically); and, most importantly for the hypothesis, **the overlap
between the model's best-projecting frequencies and the stated vocal tonic
ranges (§4.4) — if anything the alignment is slightly tighter with the
corrected cone geometry** (185Hz best point now sits inside the sean-nós
tonic range rather than at its upper edge).

This is offered to a reviewing model as the single most useful
sensitivity check already performed: a major, independently-motivated
geometry correction did not overturn the central finding.

### 4.6b Response to `review_round_notes.md`'s two "Priority 1" claims

Before acting on either, both were checked read-only against the actual
code/data (acting on a mistaken bug report would waste compute and risked
corrupting otherwise-sound results):

- **Claimed inconsistency between the headline chart and
  `ground_projection_standing.csv` at 250Hz**: not a bug. The headline
  chart correctly averages gain over all sampled azimuths, as documented;
  the reviewer's numbers were for azimuth=0 alone. This is the same
  8-azimuth undersampling issue addressed properly in §4.7 below — the
  reviewer's specific numbers were right, but the conclusion ("inconsistency")
  was a methodology mismatch, not a bug in either file.
- **Claimed hemisphere-split figure of "~49/51" vs. the reported ~35.6%
  lower**: recomputing directly from `elevation_sweep_standing.csv` with
  the same cos(elevation) solid-angle weighting already in
  `solve_tower_radiation.py` reproduces 35.6% exactly, independently
  stable to within 0.3% across all 8 sampled azimuths. A plausible
  alternative mistake (averaging dB values directly instead of converting
  to linear intensity first) lands at 47.9%, close to the reviewer's
  claimed figure. No evidence the current hemisphere-split code is wrong;
  real evidence their proposed number is reachable via a known class of
  error. Not changed. (Caveat: see open item §5.1 — the elevation sweep's
  own azimuth resolution wasn't widened this round, so its numbers carry
  the same undersampling risk §4.7 found and fixed for the ground grid.)

### 4.7 CORRECTED: frequency response at proper (36-point) azimuth resolution

Following the `review_round_notes.md` review, `AZIMUTHS_DEG` in
`solve_tower_radiation.py` was widened from 8 points (45° spacing, an exact
submultiple of the tower's 90° window symmetry) to 36 points (10° spacing),
and the full sweep re-run. The elevation sweep's azimuth resolution was
deliberately left at 8 points (a separate constant, `ELEVATION_AZIMUTHS_DEG`)
since only the ground-projection grid was flagged.

**The result changes the headline finding materially.** Gain vs. free-field
baseline at 100m, properly resolved:
- Best frequency: **125Hz, gain ≈ −5.2dB** (previously reported: 185Hz,
  ≈−1.4dB)
- Worst frequency: **50Hz, gain ≈ −9.8dB** (previously reported: 250Hz,
  ≈−12.2dB — the low-frequency end was always weak, but wasn't previously
  the single worst point)
- The whole curve now spans a much narrower range, roughly −5 to −10dB,
  with no dramatic isolated peak or trough — a gently undulating shape
  rather than the sharp two-peak/one-trough story in §4.3.

**Why the earlier numbers were wrong, concretely (250Hz case study):** at
250Hz/100m, the true azimuthal pattern (36 points) is NOT the simple
on-axis/diagonal alternation the 8-point grid suggested. It's a repeating
90°-periodic pattern with genuine internal structure: e.g. across one
90° sector, gain runs −20.7dB (0°, on a window axis) → −5.8dB (10°) →
**−2.0dB (20°, a local peak, better than the "good" diagonal points)** →
−14.0dB (30°, another dip) → −4.8dB (40°) → −4.8dB (50°) → −14.0dB (60°,
mirrors 30°) → −2.0dB (70°, mirrors 20°) → −5.8dB (80°) → −20.7dB (90°,
next window axis). The 8-point grid sampled only 0°, 45°, 90°... and
happened to land on points that made the pattern look like a clean
two-level alternation. The properly-resolved mean at 250Hz/100m is
**−8.2dB**, not −12.2dB.

**Note-overlay comparison, redone with corrected data** (same method as
§4.4, same tune profiles, not repeated here):

| Note | Freq | Gain (corrected) | Gain (original, §4.4, now superseded) |
|---|---|---|---|
| D/Do tonic | 196Hz | −7.7dB | −0.65 to −1.7dB |
| D/Do tonic | 207.65Hz | −7.2dB | −1.2 to −3.1dB |
| F minor 3rd (Martyrs) | 233.1Hz | −6.6dB | −8.1 to −10.2dB |
| F minor 3rd (Martyrs) | 246.9Hz | −8.2dB | −11.6 to −12.9dB |
| Mi major 3rd (French) | 246.9Hz | −8.2dB | −7.9 to −12.9dB |
| Sol 5th (French) | 293.7/311.1Hz | −7.4dB | −2.5 to −3.6dB |

**Honest read of this table**: the "tonic clearly beats the third" story
does not survive. At the G3 tonic assumption (196Hz), the tonic is now
actually WORSE than Martyrs' minor 3rd at 233.1Hz (−7.7dB vs. −6.6dB) —
the direction of the effect *reverses* for that specific pitch pairing.
Only under the Ab3 tonic assumption does tonic (−7.2dB) still beat the
3rd (−8.2dB), and only by 1dB, not the 5-10dB previously claimed. French's
fifth (previously the standout "second pillar aligns too" result) is now
unremarkable — solidly mid-pack, not a peak.

**What this means for the project's central claim:** the strong version of
"the tower's acoustic sweet spot lines up with the vocal tonic" is not
supported by the corrected data. A weaker, more defensible version may
still hold — the tower doesn't show strong frequency-selective behavior at
all in this range once properly resolved, so there's no sharp mismatch
either — but the striking alignment that drove several turns of discussion
was substantially a sampling artifact. This should be stated plainly to
anyone relying on the earlier framing.

### 4.8 Reciprocity check and CFIE cross-check (Tasks 5 & 6)

**Reciprocity check** (`check_reciprocity.py`): swapped source and receiver
between the chamber centre (0,0,1.75m) and a point near window 0
(1.9,0,0.7m) at 185Hz, using the validated interior rigid-cavity
formulation. Reciprocity theorem requires the transfer pressure to be the
same either way. Result: **2.6% magnitude difference, 2.3° phase
difference** — passes a 5%/5° tolerance cleanly. This is independent
evidence the interior formulation (§3.3) doesn't have a sign/jump-relation
bug, on the actual tower geometry rather than just the sphere validation
case.

**CFIE-style cross-check at 250Hz** (`solve_cfie_250hz.py`): rather than
build a full Burton-Miller CFIE (which needs the hypersingular operator —
judged too risky to introduce unvalidated purely for a diagnostic check,
see the script's docstring for the reasoning), this used an INDIRECT
single-layer formulation of the same exterior Neumann problem — a
different integral equation with different irregular frequencies than the
direct method, built only from operators (V, K') already validated
elsewhere in this project. Validated first against the closed-form
pulsating-sphere solution (0.58% error) before trusting it on the tower.

Result at 250Hz/100m:
| Azimuth | Direct method | Indirect (single-layer) | Difference |
|---|---|---|---|
| 0° (on window axis) | −20.7dB | −11.1dB | 9.6dB |
| ~40° (near-diagonal) | −4.8dB | −1.6dB | 3.3dB |

**Honest interpretation**: the two methods agree qualitatively (on-axis
substantially worse than diagonal, same direction, at both test points) —
this is the sharp directional interference structure that's genuinely
present. But they disagree by 3-10dB on the exact magnitude, despite both
independently passing clean validation (<1% error) against the closed-form
sphere test. That combination — agreement on the qualitative pattern,
disagreement on magnitude, both individually validated — is consistent with
250Hz sitting in a numerically delicate zone for at least one of the two
formulations (plausibly connected to the "irregular frequency" weakness of
the direct method flagged since §3.2, though this doesn't prove which
method is closer to the true answer). **Conclusion: treat the exact depth
of any feature near 250Hz as uncertain by at least several dB — the
qualitative existence of sharp directional nulls near window axes is
real and reproduced by two independent formulations, but neither the
original −12.2dB, the corrected −8.2dB average, nor the indirect method's
own numbers should be taken as precise.**

### 4.9 Crowd-clearing / height effect — rebuilt on basic scattering physics

§4.7's correction removed the sharp resonance/tonic-alignment story and
left a modest, fairly flat structural PENALTY (tower vs. free-field,
roughly −5 to −10dB) with no long-range amplification. That shifted the
project's centre of gravity away from "does the tower project chant over
great distances" and toward a narrower, more practical question: **does
the tower's elevation help a solo voice be heard over/through a crowd
actually gathered nearby**, e.g. during a liturgical gathering in the 2nd
ring. This section has been rebuilt twice within the same review round —
worth narrating honestly rather than presenting only the final version,
since the process itself is the useful part.

**Attempt 1 (original)**: borrowed the auditorium "seat-dip effect"
(15–20dB excess attenuation for sound grazing over SEATED audience rows,
peaking ~150Hz) as an analogy for a standing crowd — always flagged as the
weakest empirical input in the project, since seated-auditorium acoustics
is a different scattering geometry from a standing crowd.

**Attempt 2**: replaced that with a directly relevant source the user
supplied — Elena Shabalina, *The Propagation of Low Frequency Sound
through an Audience* (PhD thesis, RWTH Aachen, 2013) — which measures
STANDING crowds specifically and finds, via three independent methods
(analytical porous-medium model, BEM simulation, live festival
measurements), that sound decays LESS through a crowd than through free
air at 20–120Hz (impedance-mismatch reflection, not absorption — human-
body diffuse-field absorption at these frequencies measured at just
0.025–0.07 m²/person). That finding is real and correctly read. But
implementing it produced a nonsensical result: a confident, large ground-
singer ADVANTAGE over the tower, because Shabalina's reflection-
coefficient formula was extended across the whole vocal range and applied
only to the ground singer's path, while the tower's path — a geometry her
thesis never studied (Ch.9 explicitly flags oblique/downward incidence as
untested) — silently defaulted to zero effect. Giving one side of a
comparison an evidenced correction while leaving the other at an
unjustified default isn't a finding, it's an artifact — caught by the user
on the grounds that a ground singer beating a 30m tower through a packed
crowd contradicts basic, well-established acoustic practice (why concert
PA systems fly mid/high loudspeakers above a crowd while only bass
subwoofers stay grazing at ground level).

**Attempt 3 (current)**: rebuilt on basic, frequency-general scattering
physics instead of extrapolating a narrow low-frequency study.
`crowd_height_effect.py` now models a standing person as a scattering
obstacle of radius ~0.25m. Whether a sound wave "sees" that obstacle is
governed by the standard Rayleigh-to-geometric scattering transition for a
cylinder, controlled by the dimensionless parameter `ka` (k = wavenumber,
a = body radius): at `ka` ≪ 1 (wavelength ≫ body size — roughly
Shabalina's own 20–120Hz band) bodies are acoustically near-transparent,
matching her measured near-zero effect there; at `ka` approaching and
exceeding 1 (wavelength comparable to or smaller than body size — this
starts around 200Hz for a 0.25m body, within this project's own vocal
range) bodies become real obstacles, and a horizontal path must cross a
sequence of them. Per-person scattering cross-section is modelled as a
simple, standard smooth interpolation between the two known asymptotic
limits of cylinder-scattering theory (Rayleigh ∝ (ka)³ at small ka,
saturating to the geometric-optics limit 2a at large ka) — the same
general physics used for excess attenuation through forests/tree trunks
(cf. Price, Attenborough & Heap 1988, cited in Shabalina's own
bibliography) and treated formally by multiple-scattering theory (Twersky
1962/1966, also in her bibliography); no claim is made to those sources'
precise coefficients, only to the same well-established qualitative
physics. Crowd concentration kept at **0.5 pers/m²** (a "loose gathering"
estimate, anchored to Shabalina's own lowest live-measured festival
bracket, 0.345 pers/m²), and the in-layer path length capped at a
plausible continuous-crowd depth (20m — Shabalina's own real crowd depths,
both BEM and live, were in the 10–30m range) rather than assuming a solid
packed crowd stretches to arbitrary distances.

**Critically, the same formula is applied identically to both paths.** The
ground singer's horizontal path crosses a long sequence of bodies and
picks up the full penalty; the tower's near-vertical path crosses at most
a very short stretch of crowd near the listener and picks up very little
— by the same physical mechanism (cumulative shadowing along the direction
of travel), not by one side being measured and the other defaulted to
zero. This removes the asymmetry that broke Attempt 2. Sanity check (also
run in `crowd_height_effect.py`): at Shabalina's own validated low end
(50–80Hz) the model predicts a near-zero penalty for both paths, a genuine
cross-check against her measured result even though her specific formula
is no longer used; extended informally up toward speech-intelligibility
frequencies (500–4000Hz, well beyond this project's current vocal-range
scope), the model saturates at a stable ~20dB tower advantage — the
"tower has to win by a long margin" regime expected on ordinary acoustic
grounds once wavelength is well below body size.

**Result**: the crossover from "tower disadvantaged" (dominated by the
§4.7 structural penalty, since the crowd correction is still small) to
"tower favoured" (crowd correction grows past the structural penalty)
lands right around the middle of the practical vocal range. Within the
sean-nós tonic range specifically (131–196Hz) the two mechanisms roughly
offset — mean net effect close to zero (2nd ring: −1.2dB; 3rd ring:
−0.5dB, both averaged across azimuth and this narrow band, with real
azimuthal spread from the tower's own interference structure, roughly
−23 to +9dB). Above the tonic range (250–390Hz, still within this
project's modelled frequency grid) the tower pulls ahead more clearly,
and **more so at the 3rd ring than the 2nd** (2nd ring mean +1.6dB, 3rd
ring mean +6.8dB) — consistent with the basic-physics expectation that a
clear overhead line matters more for listeners deeper in the crowd, who
would otherwise have to look through more intervening bodies. This is a
qualitatively sensible, monotonic-with-frequency result, not the earlier
version's implausible reversal — though it should be read as illustrating
the right *shape* of the effect from first-principles scattering physics,
not as a precisely calibrated number (the cross-section interpolation,
concentration estimate, and depth cap are all simple, order-of-magnitude
choices, not fitted or independently measured for this exact scenario).

**What this doesn't cover**: this project's practical vocal range has so
far been capped around 110–390Hz (the sean-nós fundamental/tonic range).
Speech *intelligibility* proper depends heavily on consonant/sibilant
energy well above that (roughly 1–4kHz), where this model's own numbers
suggest the tower's advantage would be substantial and clear-cut — but
extending the project's modelled frequency range that far, and revisiting
what "intelligibility" would even mean for chant rather than speech, is
future work, not attempted here.

**Directional structure, not a uniform circle**: re-running the advantage
calculation WITHOUT averaging over azimuth (`plot_advantage_zone_directional.py`)
still shows the same 4-fold/8-lobe interference structure found in
§4.7/§4.8, now layered on top of a crowd correction that shifts from
mildly negative to mildly positive depending on frequency and distance
rather than uniformly favouring either side — at 165Hz specifically the
map is a genuine mix of green and red, not dominated by either.

**Emerging hypothesis, explicitly unverified, and unaffected by the crowd-
model revision above** (since it concerns the tower's own structural
interference lobes, not the crowd correction): does this directional
structure show up in real round tower enclosure shapes/positioning — e.g.
enclosure boundaries or gathering-area concentration preferentially
following the cardinal/intercardinal directions from the tower, rather
than a simple one-sided elongation? A literature search turned up general
confirmation of the standard site-layout convention this project already
assumed (tower west of church, door facing the church doorway) and
confirmation that enclosures do vary between circular and oval (40–400m),
but **no existing study testing any correlation between enclosure shape
and tower acoustics** — this appears to be a genuinely novel angle. See
§4.10 for a first, partial, structural-only (not crowd-corrected) check of
a narrower version of this question against real site data. Testing the
broader enclosure-*outline* version needs comparative data this project
does not have — flagged as an open item (§5), not pursued further here.

Diagrams (all regenerated for the basic-physics model): `height_mechanism_diagram.png`
(ray-geometry explainer — why the ground singer's path crosses a long
sequence of bodies and the tower's doesn't), `advantage_zone_map.png`
(azimuth-averaged), `advantage_zone_directional.png` (directionally-
resolved — use this one, not the averaged one, for the real shape of the
advantage/disadvantage zone, which now depends on frequency and no longer
uniformly favours either side), `crowd_height_effect.png` (frequency
sweep at pilgrim-ring distances, showing the crossover from tower-
disadvantaged to tower-favoured as frequency rises).

### 4.10 Real comparative site data (Herity 1984) — the enclosure-shape hypothesis gets sharper

The user supplied a real, citable source: Michael Herity, "The layout of
Irish early Christian monasteries," in *Ireland and Europe: The Early
Church* (Próinséas Ní Chatháin & Michael Richter, eds., Stuttgart:
Klett-Cotta, 1984), pp. 105–116. It gives quantified tower-to-church
positions at three round-tower sites:

| Site | Tower-to-church distance | Herity's stated direction | Geometrically-implied bearing | Notes |
|---|---|---|---|---|
| Clonmacnois | 30m | NW of west doorway | NW (stated directly) | "its doorway facing that of the church" |
| Kells (Ceannanas Mór) | 25m | "S" of church | **SW** (corrected, see below) | doorway "facing across the open space to the west façade" |
| Glendalough | 15m | NW of cathedral west façade | NW (stated directly) | doorway "facing across the open space to the cathedral" |
| Slane | (unspecified) | "in the arc west" of church | west-ish, unspecified | speculative reconstruction, not firmly established |

**Correction to the Kells reading**, from the user, geometrically sound:
Herity's text says the tower is "to the south" of the church AND that its
doorway faces "across the open space to the *west façade*" specifically.
A tower due south of the church would face the church's south wall, not
its west façade — for the doorway to target the west façade across the
shared open space (*platea*), the tower must be offset **south-west**,
not due south. So the corrected reading across all three quantified
sites: the general rule of westward siting is adhered to, but the precise
bearing is consistently offset (NW at Clonmacnois/Glendalough, SW at
Kells) rather than exactly due west, while the tower door still targets
the church across the shared open space. This directly confirms §2's
caveat that this project's own "exactly due west" convention was a
simplification, not a documented general rule.

**A sharper, checkable hypothesis follows from this.** If the four
windows stay cardinal (N/E/S/W) regardless of where the tower sits
relative to the church — which the surviving evidence for window
orientation generally supports, though this project has no site-specific
confirmation either way — then an NW or SW tower-to-church offset means
that bearing is **intercardinal** (NE/SE/SW/NW): it falls exactly between
two windows, not in line with either one. That's directly checkable
against this project's own data (`ground_projection_standing.csv`, all 34
frequencies, 36 azimuths — no new solve needed, purely a re-read of
existing results):

- Averaged across all 34 frequencies, diagonal (intercardinal) azimuths
  beat on-axis (window-facing) azimuths overall: **−4.55dB vs. −5.78dB**.
- **Not a uniform win.** Diagonal decisively rescues the worst bands
  (+10 to +14dB better right around the 220–250Hz trough; +7 to +16dB
  better around 75–100Hz), but on-axis is actually 3–4dB *better* in the
  140–190Hz band — which is the model's best-performing range overall.
- Net reading: an intercardinal tower-to-gathering-space bearing doesn't
  universally help. It specifically buys insurance against the worst
  frequency-dependent outcomes (including the trough), at a modest cost
  in the already-good band. Since the practical vocal range spans both
  regimes, this leans toward intercardinal siting being the safer choice
  on balance — not a clean universal win, and this should be stated with
  that precision, not oversold as a tidy confirmation.

**What this does and doesn't establish**: this is a real, specific,
internally-consistent pattern — three real sites, all intercardinally
offset in a way that (if the cardinal-window assumption holds) lands the
church/gathering direction in a mostly-favourable acoustic zone relative
to this project's own interference model. But it rests on n=3 site
positions and an unverified assumption about window orientation at those
specific towers, and doesn't touch the broader enclosure-*shape* version
of the hypothesis (§4.9), which still needs outline-level data this
project doesn't have (deferred, along with a search for the original
survey sources Herity cites — Henry 1957/65/67/70, Leask 1955, Macalister
1909, Wakeman 1893 — some of which may be public-domain and more
detailed than Herity's own summary schematics). Treat this section as a
sharpened, partially-supported hypothesis, not a finding.

**Reconstructed core layouts, including crosses and other monuments**
(`plot_reconstructed_site_plans.py`, `reconstructed_site_plans.png`):
Herity's text gives distance+bearing from the church for several other
monuments at these same three sites, letting each site's small cluster of
church/tower/cross(es) be reconstructed approximately (positions computed
from the stated vectors, not digitized from a survey). Cross-to-tower
distances are **derived by triangulation** — Herity states cross-to-church
and tower-to-church vectors separately, never cross-to-tower directly:

| Site | Monument | Distance + bearing from church (stated) | Distance from tower (derived) |
|---|---|---|---|
| Clonmacnois | South Cross | 12m S, 4m E of façade | 41.7m |
| Clonmacnois | North Cross | 12m N, 4m E of façade | 26.8m |
| Clonmacnois | Cross of the Scriptures | 12.5m from west façade (direction assumed W) | 22.9m |
| Clonmacnois | Teampall Ciarán | NE of church (no distance stated) | — |
| Kells | Cross of Patrick & Columba | "beside" the tower — not independently measured | ~0m (co-located) |
| Kells | Broken Cross | opposite the west doorway (no distance stated) | — |
| Kells | St Columb's House | 90m NW of west façade | 93.4m |
| Glendalough | plain granite High Cross | south of façade (no distance stated) | — |
| Glendalough | Priest's House | 27m SSW, on a lower terrace | 35.6m |

At Clonmacnois specifically, all three crosses sit farther from the
church (12.5–12.7m) than they do from the tower's own distance scale
(22.9–41.7m from the tower vs. the tower's 30m from the church) — the
crosses cluster tightly around the church, not around the tower or
between the two. That's a useful check against over-reading the
tower-bearing pattern in §4.10 above: the crosses' positions look driven
by proximity to the church/open-space focus described throughout
Herity's paper, not by any relationship to the tower specifically. Worth
being explicit that this doesn't strengthen or weaken the diagonal-bearing
acoustic hypothesis either way — it's a different, independent siting
logic (crosses cluster on the church's open space; the tower sits at its
edge) that this data happens to let us see clearly for the first time.

## 5. Open items / recommended checks for a reviewing model

**Resolved this round** (see §4.7/§4.8): the 8-azimuth aliasing risk was
real and materially changed the headline result; reciprocity check passed;
CFIE-style cross-check gave a nuanced (qualitative-agreement, quantitative-
disagreement) result at 250Hz. Two follow-up items explicitly deferred by
the user given compute cost: TOWER_RADIUS sensitivity sweep (1.9m/2.5m,
~100 min), and a window radiation-impedance coupling check (lowest priority
per the review document itself).

1. **The elevation sweep still uses 8 azimuth points** (`ELEVATION_AZIMUTHS_DEG`,
   deliberately left unwidened since only the ground-projection grid was
   flagged) — given §4.7 showed real fine angular structure the 8-point
   grid missed for the *ground* grid, the elevation/vertical-directivity
   conclusions in §4.5 (the "~1/3 lower hemisphere" figure) may have the
   same undersampling risk. Worth widening and re-checking before trusting
   §4.5's specific number, even though the hemisphere-split *methodology*
   itself was independently confirmed correct (§4 Task 2 finding above).
2. **Robustness of the (now much flatter) frequency response to the
   remaining guessed parameters** (damping factor, source calibration, and
   the TOWER_RADIUS mismatch flagged in §2). Given how much the properly-
   resolved azimuth grid already changed the picture, a sensitivity sweep
   over these remaining parameters is more important now, not less —
   deferred this round, but higher priority than before this correction.
3. **The weak interior/exterior coupling approximation** (§3.3) — is
   treating ~5% window area as negligible for the interior solve
   defensible, or should a fully self-consistent transmission problem be
   attempted? Still not done (Task 7, explicitly deferred).
4. **Re-examine the "did the tower influence pitch choice" speculation**
   given §4.7's correction — the alignment it was speculating about no
   longer clearly exists. This line of reasoning should probably be
   dropped or substantially walked back, not just re-tested.
5. Sanity-check the overall claim structure now that §4.4's note-overlay
   alignment is superseded (§0) and §4.9's crowd-clearing case has been
   rebuilt on basic scattering physics rather than the original seat-dip
   analogy. The result is more mixed than either prior version claimed:
   within the sean-nós tonic range (131–196Hz) the crowd mechanism and the
   §4.7 structural penalty roughly cancel (net close to zero at both real
   pilgrim rings); above the tonic range, within the rest of this project's
   modelled frequency grid (250–390Hz), the crowd mechanism increasingly
   favours the tower, more so at the 3rd ring than the 2nd. What remains
   standing throughout: interior symmetry from source-centring (§4.1,
   robust, unaffected by any correction so far), and a general (not
   sharply frequency-selective) projection efficiency of −5 to −10dB
   relative to free-field (§4.2/§4.7 — a real, modest structural penalty).
   The project's centre of gravity has shifted from a single clean
   narrative ("the tower helps") toward a frequency-dependent picture: this
   project's own numbers now suggest height's practical value is weak to
   neutral for the low/tonic end of the vocal range and grows toward being
   a clear advantage further up — with speech-intelligibility frequencies
   (well above this project's current 390Hz ceiling) plausibly the range
   where the original hypothesis is best supported, not the tonic range
   the project has focused on so far. That's a genuinely open, useful
   question for a reviewing model to weigh in on.
6. **Crowd-clearing model (§4.9) — went through two revisions this round,
   both worth knowing about.** First, a directly relevant source (Shabalina
   2013, standing-crowd low-frequency acoustics) was found and applied in
   place of the seated-audience seat-dip figure — correctly read, but
   mis-applied: her reflection-coefficient formula was extended across the
   whole vocal range and applied only to the ground singer's path, while
   the tower's untested geometry silently defaulted to zero effect,
   producing an implausible, large ground-singer advantage that the user
   caught as inconsistent with basic acoustic practice. That version was
   discarded. The current version instead uses basic, frequency-general
   cylinder-scattering physics (Rayleigh-to-geometric transition vs. body
   size), applied identically to both paths — a smaller, more defensible
   effect that reproduces Shabalina's own near-zero low-frequency finding
   as a consequence rather than an assumption, and grows into a real
   tower-favouring effect at higher frequency. This is a simple,
   order-of-magnitude model (a standard textbook interpolation between two
   known scattering limits, not a fitted or independently validated
   number for this exact scenario) — a tighter calibration, and extending
   the frequency range further toward intelligibility-relevant frequencies,
   are candidates for a future round, not attempted here.
7. **The enclosure-shape/directional-skew hypothesis is partially
   addressed (§4.10), still open at the outline level.** Real bearing data
   for 3 sites (Herity 1984) now supports the narrower claim — towers
   offset NW/SW of the church, landing the gathering direction on the
   acoustically-favourable-on-balance intercardinal bearing in this
   project's own model — but n=3, and the window-orientation assumption
   feeding that claim is unverified at those specific sites. The broader
   enclosure-*outline*-shape version (does the enclosure boundary itself
   skew toward the good directions, not just the tower-church bearing)
   is still untested; needs enclosure outline data this project doesn't
   have. Two possible next sources, both deferred so far: the original
   survey publications Herity cites (Henry 1957/65/67/70, Leask 1955,
   Macalister 1909, Wakeman 1893 — some pre-1923 and possibly public
   domain), or digitizing Herity's own low-resolution schematic figures
   (Fig 3B/Kells and Fig 2C/Clonmacnois show outline-relevant detail,
   Fig 2D/Glendalough and Fig 3A/Armagh less so) — the latter would only
   ever be illustrative given the reproduction quality, not precise.

## 6. Software / reproducibility notes
Python + Gmsh (OpenCASCADE) for geometry/meshing, `bempp-cl` 0.4.2 for BEM
(CPU/Numba backend, no OpenCL/FMM available), NumPy/SciPy/Matplotlib.
Environment: a project-local venv (`Round_Tower_project/.venv`,
`requirements.txt` in the repository root). Key scripts: `build_tower_mesh.py`,
`validate_bem_formulation.py`, `validate_interior_bem.py`,
`solve_interior_field.py`, `solve_tower_radiation.py` (this now also
produces the pilgrim-ring-walk CSV directly, evaluated from the same
solved field as the rest of its output — `evaluate_pilgrim_rings.py` was
the original standalone version, which wastefully re-solved the whole BEM
problem a second time purely to evaluate different points; superseded,
kept only for reference), `check_reciprocity.py`, `solve_cfie_250hz.py`
(§4.8), `crowd_height_effect.py` (§4.9, post-processing only, no BEM),
plus various `plot_*.py` scripts including `plot_height_mechanism.py` and
`plot_advantage_zone_directional.py` (§4.9 diagrams).
