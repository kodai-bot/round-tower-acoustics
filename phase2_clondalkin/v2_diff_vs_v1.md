# v2 vs v1 headline diff

**All v1 frequency-domain numbers below are SUPERSEDED.** v1 modelled the chamber as a zero-thickness shell at a single compromise radius (2.2m, a general cross-tower average) -- roughly 2x the true Clondalkin internal diameter, with wall thickness omitted entirely and windows as flush patches rather than real ducts. v2 rebuilds on RIA-drawing-derived, Clondalkin-specific dimensions with a real ~0.745m wall and rectangular duct windows. Nothing below should be read as a correction to v1's numbers -- it is a different, more physically-grounded model of a specific real tower.

## Headline projection efficiency (gain vs. free-field @ 100m)

| | v1 (superseded, zero-thickness shell) | v2 (Clondalkin, real duct wall) |
|---|---|---|
| Best frequency | 125Hz, -5.2dB | 75Hz, 17.7dB |
| Worst frequency | 50Hz, -9.8dB | 125Hz, -13.5dB (converged frequencies only) |

**Read the two caveats below before trusting either number in this table.**

### Caveat 1: numerical reliability boundary -- not every point in the 45-point sweep is trustworthy

GMRES was capped at 400 iterations (`v2_params.GMRES_MAXITER`, added after an
uncapped first attempt at this sweep ran 10+ hours -- see the sweep log for
the full story). 18 of 45 frequencies hit that cap without converging to the
1e-5 tolerance: **250, 260, 280, 290, 300, 310, 320, 340, 349, 360, 375, 390,
400, 420, 440, 460, 480, 500 Hz**. The clean, converged set is the other 27:
**50, 59, 75, 100, 118, 125, 143, 150, 160, 165, 170, 175, 180, 185, 190, 196,
200, 207, 210, 220, 230, 237, 240, 270, 293, 296, 330 Hz**.

The pattern is not a clean "above 300Hz" cutoff -- 270/293/296/330Hz all
converged cleanly, scattered among the failures -- but the failure rate is
overwhelming from 300Hz up (14 of 15 points in that range failed, only
330Hz converged) and already substantial in the 250-290Hz decade (4 of 5
points failed). The worst-case figure in the table above is restricted to the 27 clean
points: 125Hz, -13.5dB. (An earlier version of this table reported 440Hz,
-15.0dB, which came from a non-converged solve; `v2_diff_vs_v1.py` now
excludes `v2_params.NONCONVERGED_HZ`.)
This mirrors, but is more severe than, v1's own long-flagged caveat that
the direct (non-CFIE) BIE formulation can fail near "irregular
frequencies" coinciding with interior Dirichlet eigenvalues -- here it
looks like the mesh-resolution margin (0.15m elements, chosen as a
compute-budget fallback -- see `v2_params.MESH_SIZE`'s comment) compounds
with that known weakness as frequency rises, rather than being purely an
isolated-frequency effect.

### Caveat 2: the +17.7dB "best frequency" at 75Hz is real but its exact size is uncertain

75Hz DID converge cleanly (not one of the 18 flagged above), and its result
is numerically stable -- re-solved at 1000x tighter GMRES tolerance
(1e-8 vs 1e-5), the answer changed by less than 0.001% (`|p_tower|`
2.835464e-3 vs 2.835462e-3). So this is not a GMRES artifact.

It was still cross-checked against the independent indirect (CFIE-style)
formulation used elsewhere in this project (see `v2_solve_cfie_crosscheck.py`),
because a sharp, narrow low-frequency peak in an *exterior* radiation
problem is an unusual enough shape to warrant the same scrutiny v1 gave
its own sharp features. **Both methods agree there is a real, substantial
positive gain at 75Hz** -- the direct method gives +17.75dB, the indirect
method +9.41dB (at 0deg/100m) -- qualitatively the same striking result
(nothing else in the neighbouring frequencies comes close to positive
gain at all), but an ~8dB quantitative gap between the two independently-
formulated methods. Report this as "a genuine, cross-validated low-
frequency enhancement effect in the roughly +9 to +18dB range," not as a
precise +17.7dB figure -- consistent with this project's established
practice (v1's 250Hz CFIE cross-check) of reporting qualitative agreement/
quantitative uncertainty honestly rather than picking one number.

The 190Hz point (this sweep's single largest azimuth-to-azimuth swing,
23.9dB at 100m) was also cross-checked: direct and indirect again agree
on which direction is favoured (on-axis better than the 45deg diagonal)
but disagree substantially on magnitude at the diagonal point (-22.91dB
direct vs -2.09dB indirect, a 20.8dB gap) -- the same qualitative-yes/
quantitative-uncertain pattern, in a different part of the spectrum.

## Four-window balance (centred standing source)

| | v1 | v2 |
|---|---|---|
| Mean imbalance across sweep | 0.02dB | 0.09dB |
| Worst-case imbalance | 0.17dB | 1.75dB |

## Vertical directivity (power-weighted lower-hemisphere fraction)

| | v1 | v2 |
|---|---|---|
| Fraction radiated below horizontal | 35.6% | 43.6% |

(v1 elevation grid: 8 azimuths, unchanged this round in v1 -- see PROJECT_SUMMARY_FOR_REVIEW.md SS5 item 1. v2 elevation grid: 16 azimuths.)

## What changed in the geometry itself

| Parameter | v1 | v2 |
|---|---|---|
| Interior radius | 2.2m (general average) | 0.9905m (Clondalkin RIA drawing) |
| Wall thickness | 0 (zero-thickness shell) | 0.745m (real, ducted) |
| Window model | flush Dirichlet patch | rectangular duct through full wall depth |
| Chamber height | 1.2192m (cylinder) + 2.5908m cone | 2.896m (Variant A, ceiled -- cone excluded from acoustic domain) |
| Source position | (0,0,1.75m), above window band | (0.0, 0.0, 1.5), inside window band |
| Frequency ceiling | 390Hz | 500Hz |
| Ground azimuth resolution | 36 points (widened mid-round after an aliasing bug) | 36 points (from the start) |
| Mesh topology | open tube, missing floor cap (found during v2 planning, not corrected in v1 -- see v2_check_watertight.py) | closed, watertight, checked automatically |
