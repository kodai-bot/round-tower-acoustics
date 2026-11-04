# Round Tower BEM pipeline — next round of checks

Context: this is a follow-up to the "more robustly constructed model" summary (cone-roof
geometry, validated formulations). A reviewing pass over the results surfaced two
apparent inconsistencies that should be reconciled before the numbers are trusted
further, one known undersampling issue, and several robustness checks that were
already flagged as open items. Tasks are ordered by priority; each lists what to do,
why, and what output to produce.

---

## Priority 1 — Reconcile two inconsistencies

### Task 1: Confirm what the headline frequency-response sweep is actually evaluating
- Find the script/function that produced `headline_frequency_response_cone_labeled.png`
  (likely inside `solve_tower_radiation.py`).
- Determine whether it evaluates pressure at a literal fixed point (r=100m, azimuth
  0°/east, via the potential operator) or via bempp's asymptotic far-field operator
  (`far_field.helmholtz.single_layer`, as used elsewhere in the codebase).
- Cross-check against `ground_projection_standing.csv` at azimuth=0: the 250Hz gain
  there ranges from about −14dB (50m) to −20.7dB (100m) to −15.1dB (200m), settling
  to roughly −13.7 to −13.9dB beyond 500m. The headline table reports −12.2dB at
  250Hz. Work out which evaluation method + distance combination actually reproduces
  −12.2dB (the converged large-distance value is the closest match found so far).
- Deliverable: a short note (comment block or README) stating precisely what
  "gain @ 100m" in the headline plot means operationally, and why it agrees closely
  with the 185Hz point-evaluation but not the 250Hz one (this may itself be a
  meaningful finding — 250Hz sits near a sharp, spatially unstable interference null,
  185Hz sits on a broad stable peak).

### Task 2: Reconcile the elevation upper/lower hemisphere split
- Locate whatever produced the "~1/3 lower hemisphere / 2/3 upper" figure quoted in
  the review summary (§4.5) — likely derived from `farfield_sphere.vtu` /
  `farfield_pattern_db.txt` or an earlier ad hoc calculation.
- Compare its methodology (solid-angle weighting, angular resolution, evaluation
  distance) against `elevation_sweep_standing.csv`, which — properly solid-angle
  weighted with cos(elevation) — gives a near-even ~49% lower / 51% upper split,
  consistently across all 8 sampled azimuths and across frequency.
- Work out which figure is more reliable (likely whichever has finer angular
  resolution and correct solid-angle weighting) and update the write-up to match, or
  explain clearly why both are legitimate under different definitions.

---

## Priority 2 — Fix a known undersampling risk

### Task 3: Re-run the ground-projection sweep at finer azimuth resolution
- Same distance grid (50, 100, 200, 500, 1000, 2000m) and same 34-point frequency
  sweep as the current `ground_projection_standing.csv`.
- Increase azimuth sampling from 8 points (45° spacing) to at least 36 points (10°
  spacing), matching the resolution already used in `pilgrim_ring_walk_standing.csv`.
  Reason: the tower has exact 4-fold (90°) symmetry, so 45°-spaced sampling is a
  submultiple of that symmetry and risks systematically aliasing with the real lobe
  structure rather than randomly missing it — confirmed by comparing this file's
  185Hz result (looks nearly omnidirectional beyond 100m) against the finer
  pilgrim-ring data at the same frequency (shows ~20dB of real angular structure).
- Output: `ground_projection_standing_fine.csv`, same column structure as the current
  file.

---

## Priority 3 — Robustness checks (carried over from the original review)

### Task 4: TOWER_RADIUS sensitivity sweep
- Re-run the full mesh-generation + solve pipeline at chamber radius = 1.9m and 2.5m
  (in addition to the current 2.2m), holding roof/window geometry (from the
  Clondalkin scan) fixed.
- Record how the headline curve's peak (~185Hz) and trough (~250Hz) frequencies and
  gains shift with radius.
- Output: small CSV/table of `{radius, peak_freq, peak_gain, trough_freq, trough_gain}`.

### Task 5: Reciprocity check on the two boundary-integral formulations
- At one test frequency (185Hz is fine): solve once with a monopole source at a
  window location and evaluate pressure at the chamber centre; solve again with
  source and receiver swapped. The two transfer pressures should agree (reciprocity
  theorem) — this is a cheap independent check for sign/jump-relation errors in
  either formulation, without a full re-derivation.
- Output: printed comparison + pass/fail note.

### Task 6: Direct vs. CFIE cross-check specifically at 250Hz
- At 250Hz only, solve the exterior problem with the current direct (non-CFIE)
  formulation and also with a CFIE/Burton–Miller combined formulation (alpha ≈ 0.5,
  same form as the original `bempp_radiation.py` scaffold), to test whether the
  trough is affected by proximity to an irregular frequency of the direct method.
- Output: gain_db from both formulations at 250Hz, church-facing bearing, a few
  distances.

### Task 7 (lower priority): Window radiation-impedance coupling check
- At 1–2 representative frequencies, replace the fully-rigid window boundary
  condition in the interior solve with a simple piston-in-infinite-baffle radiation
  impedance (function of ka), and compare the resulting window-plane pressure
  against the current rigid-wall approximation. Quantifies how much the "5% window
  area is negligible" simplification (§3.3 of the review) actually affects the
  driving pressure.

---

## Notes for the agent
- Tasks 1–3, 5–7 reuse existing mesh/solve infrastructure; only Task 4 needs new
  geometry variants.
- Please keep output file naming consistent with the existing convention
  (`..._standing.csv`, `..._labeled.png`) and note in each output which script
  version produced it, given this project's history of superseded scripts.
