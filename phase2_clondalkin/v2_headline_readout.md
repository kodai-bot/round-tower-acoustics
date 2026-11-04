# Round Tower v2 (Clondalkin) — Headline Results Readout

Purpose of this file: a self-contained reference for discussing "what would
this sound like at distance X" and other headline v2 findings with another
agent/tool, without needing to re-derive numbers from the raw CSVs. Every
number below is sourced from a specific v2 output file (paths given) so it
can be checked or extended.

**Read the "Reliability" section before using any number below.** Not every
frequency in this sweep converged cleanly — using a flagged-unreliable
point without noting that would misrepresent the model.

---

## 0. What this model is (one paragraph)

A Boundary Element Method (BEM) acoustic model of the upper chamber of an
Irish Round Tower, rebuilt this round on measurements specific to
**Clondalkin** (Ireland's narrowest surviving round tower), from an RIA
archive cross-section drawing: real ~0.745m-thick stone wall, 4 cardinal
windows modelled as actual rectangular ducts through that wall (not flush
holes), rigid ceiling sealing the chamber below the roof cone (the
"Variant A" / ceiled case). Source: a single standing voice on the chamber
axis, calibrated to 84dB SPL at 1m (a "trained singer" level assumption,
unchanged from round 1 — no measurement exists for real sean-nós/precentor
delivery). Tower floor assumed 30m above ground. All SPL figures below
follow directly from that 84dB@1m calibration — if you want to reason
about a louder or quieter voice, every SPL number scales linearly in dB
with the source level (e.g. a 90dB@1m voice → add 6dB to every SPL figure
below; gain_db figures are unaffected, they're ratios).

Full method/geometry writeup: `rebuild_spec.md` (spec) and
the docstrings of `build_tower_mesh_v2.py`, `solve_interior_field_v2.py`,
`solve_tower_radiation_v2.py`. Full diff against the superseded round-1
model: `v2_diff_vs_v1.md`.

---

## 1. Reliability — which frequencies to trust

GMRES was capped at 400 iterations per solve. **27 of the 45 swept
frequencies converged cleanly; 18 did not** (hit the cap without reaching
1e-5 tolerance) and should be treated as unreliable, not as real findings:

- **Reliable (27):** 50, 59, 75, 100, 118, 125, 143, 150, 160, 165, 170,
  175, 180, 185, 190, 196, 200, 207, 210, 220, 230, 237, 240, 270, 293,
  296, 330 Hz
- **Unreliable / non-converged (18), do not quote these:** 250, 260, 280,
  290, 300, 310, 320, 340, 349, 360, 375, 390, 400, 420, 440, 460, 480,
  500 Hz

Practical upshot: **the model's effective trustworthy range is ~50–330Hz**,
with a weak spot at 250–290Hz. This comfortably covers the sean-nós vocal
tonic range this project cares about (roughly 110–250Hz, sweet spot
130–196Hz per round-1's user-supplied vocal-range data), but the 300–500Hz
extension the rebuild instruction asked for did not end up delivering
trustworthy data — that's a real, load-bearing limitation, not a footnote.

Everything in the tables below uses only reliable frequencies.

---

## 2. SPL at distance, by frequency — the "what would it sound like" table

Absolute SPL (dB), azimuth-averaged, at the given straight-line ground
distance from the tower base. Source: `v2_ground_projection_standing.csv`.
For reference: ~20dB ≈ quiet rural night, ~40dB ≈ quiet room/library,
~50-60dB ≈ normal conversation at 1m, ~85dB starts to be genuinely loud
close-up. These are tonal single-frequency SPLs (a chant's fundamental
only, not full timbral loudness perception).

| Freq | 15m | 25m | 30m | 50m | 100m | 200m | 500m | 1000m | 2000m |
|---|---|---|---|---|---|---|---|---|---|
| 75 Hz *(see §4)* | 46.9 | 46.9 | 46.6 | 44.8 | 39.9 | 33.9 | 25.7 | 19.6 | 13.5 |
| 125 Hz (low tonic) | 58.5 | 55.3 | 53.4 | 45.4 | 34.1 | 26.3 | 18.7 | 13.2 | 7.5 |
| 150 Hz (sweet spot) | 68.5 | 64.8 | 62.6 | 49.0 | 46.8 | 39.8 | 30.0 | 24.4 | 18.9 |
| 165 Hz (sweet spot) | 47.3 | 42.0 | 40.3 | 35.7 | 34.6 | 27.1 | 18.8 | 13.3 | 7.6 |
| 185 Hz | 64.1 | 52.1 | 45.8 | 41.5 | 44.8 | 45.2 | 38.7 | 32.8 | 26.8 |
| 196 Hz (high tonic) | 48.9 | 54.6 | 53.1 | 48.5 | 45.8 | 39.8 | 32.1 | 26.3 | 20.3 |

**Note the large swings even at fixed distance as frequency changes (e.g.
150Hz vs 165Hz at 15m: 68.5 vs 47.3dB, over 20dB apart)** — this is the
tower's real, confirmed 4-fold directional interference structure (see
round-1's §4.7/§4.8) landing differently at each frequency, not noise.
Position around the tower matters as much as distance does.

15/25/30m are included because they're the real tower-to-church distances
at three excavated sites (Clonmacnois/Kells/Glendalough, round-1's Herity
1984 addendum) — useful for "how would this have sounded at the actual
gathering point" questions.

### Same table as gain vs. free-field (dB) — how much the tower structure itself helps/hurts vs. an unobstructed source of equal power

| Freq | 15m | 25m | 30m | 50m | 100m | 200m | 500m | 1000m | 2000m |
|---|---|---|---|---|---|---|---|---|---|
| 75 Hz | +14.5 | +15.9 | +16.4 | +17.4 | +17.7 | +17.4 | +17.1 | +17.0 | +16.9 |
| 125 Hz | +0.7 | -1.1 | -2.3 | -7.4 | -13.5 | -15.5 | -15.3 | -14.8 | -14.5 |
| 150 Hz | +16.3 | +14.1 | +12.6 | +1.9 | +4.8 | +3.6 | +1.7 | +2.1 | +2.6 |
| 165 Hz | -5.5 | -9.4 | -10.3 | -12.1 | -8.1 | -9.8 | -10.2 | -9.7 | -9.3 |
| 185 Hz | +5.2 | -5.3 | -10.8 | -12.2 | -3.9 | +2.3 | +3.7 | +3.8 | +3.8 |
| 196 Hz | -13.1 | -6.0 | -6.7 | -8.4 | -6.0 | -6.2 | -6.0 | -5.8 | -5.8 |

Interesting reversal from round 1: **125Hz was round 1's single
best-projecting frequency; in v2 it's the worst of the reliable set**
(-13.5 to -15.5dB at 100m+). This is a real consequence of the geometry
rebuild (much smaller, thick-walled, ducted chamber), not a bug — flagged
here because it's the kind of "don't assume round 1's specific-frequency
findings still hold" trap worth being explicit about.

---

## 3. At the real pilgrim gathering rings

Mean and [min, max] SPL across all 36 sampled positions around each ring.
2nd ring = 45m radius (accounting for the tower's 17m off-centre site
position), 3rd ring = 91.5m radius. Source: `v2_pilgrim_ring_walk_standing.csv`.

| Freq | 2nd ring mean [range] | 3rd ring mean [range] |
|---|---|---|
| 75 Hz | 45.1 [43.6, 46.8] dB | 40.6 [39.3, 42.3] dB |
| 125 Hz | 47.0 [41.4, 54.2] dB | 35.4 [32.9, 38.3] dB |
| 150 Hz | 51.9 [41.8, 63.6] dB | 46.7 [45.8, 47.7] dB |
| 165 Hz | 37.3 [33.4, 41.4] dB | 35.1 [32.6, 38.1] dB |
| 185 Hz | 43.6 [27.5, 54.1] dB | 44.2 [40.8, 49.7] dB |
| 196 Hz | 49.4 [41.4, 58.6] dB | 46.4 [38.8, 54.6] dB |

The wide within-ring ranges (e.g. 185Hz 2nd ring: 27.5 to 54.1dB, a 26.6dB
spread just from WHERE you stand on the ring) is the same directional
structure as §2 — worth emphasizing to anyone reasoning about "how loud
was it for the congregation": position relative to the 4 windows likely
mattered more than which ring you stood on.

---

## 4. The 75Hz anomaly — genuine but uncertain in magnitude

75Hz shows a striking, isolated **+17.7dB gain vs. free-field**, unlike
anything at neighbouring frequencies. This was investigated, not taken at
face value:
- **Numerically stable**: re-solved at 1000x tighter GMRES tolerance
  (1e-8 vs 1e-5), the result changed by <0.001% — not a solver artifact.
- **Cross-validated** against an independently-formulated BEM method
  (indirect single-layer, same technique round 1 used at 250Hz): that
  method also shows a large positive gain, +9.4dB (vs the direct method's
  +17.7dB) at 0°/100m. Both methods agree the effect is real and
  substantial; they disagree by ~8dB on its exact size.
- **Conclusion: treat this as "a genuine, cross-validated low-frequency
  enhancement effect, magnitude somewhere in the +9 to +18dB range,"** not
  as a precise +17.7dB figure. This is a genuinely new finding with no
  round-1 analogue (round 1's mesh/geometry never showed anything like it).

A second sharp feature — the largest azimuth-to-azimuth swing in the whole
sweep, 23.9dB at 190Hz/100m — was cross-checked the same way: both methods
agree on-axis beats the diagonal there, but disagree by ~21dB on the
diagonal point's exact depth. Same lesson: the qualitative directional
structure is trustworthy, exact depths in these numerically delicate
spots are not.

---

## 5. Crowd-clearing / height-advantage (tower vs. a ground-level singer)

Same cylinder-scattering model round 1 used (bodies as scattering
obstacles, Rayleigh-to-geometric transition with frequency), reapplied to
v2's geometry — see `v2_crowd_height_effect.py`. "Tower advantage" = how
many dB better the tower-mounted voice does vs. an equal-power singer at
ground level, once a crowd's scattering effect on each path is included.
Mean across all reliable frequencies and azimuths, by distance
(`v2_crowd_height_effect.csv`):

| Distance | Tower advantage (dB) |
|---|---|
| 15m | +6.7 |
| 25m | +8.6 |
| 30m | +8.5 |
| 50m | +8.1 |
| 100m | +8.8 |
| 200m | +9.7 |
| 500m | +8.5 |
| 1000m | +6.5 |
| 2000m | +2.6 |

The tower consistently helps by roughly 6-10dB across the whole practical
distance range in this model — a clearer, more one-sided result than
round 1's own finding, where the tower's advantage roughly cancelled the
structural penalty specifically within the sean-nós tonic range. (Caveat:
this table averages over reliable frequencies only, weighted equally —
it is not re-weighted by how much time a real chant spends on each pitch.)

---

## 6. Other headline single-number results

| Finding | v2 | v1 (superseded) |
|---|---|---|
| Four-window balance, mean imbalance (reliable freqs) | 0.037dB | 0.02dB |
| Four-window balance, worst-case (reliable freqs) | 0.19dB | 0.17dB |
| Power-weighted fraction radiated into lower hemisphere | 43.6% | 35.6% |
| Interior validation gate | 10/10 analytic modes matched | (not performed in v1) |
| Reciprocity check | 0.005% mag / 0.01° phase (excellent) | 2.6% / 2.3° (passed) |

The four-window balance result is essentially unchanged and re-confirms
round-1's key structural finding: a centred, on-axis standing source
drives all 4 windows almost perfectly evenly — this falls out of the
room's own near-axisymmetry, independent of the major geometry rebuild.

---

## 7. Where to look for more

- `v2_ground_projection_standing.csv` — full grid: 27 reliable freqs (+18
  flagged) × 36 azimuths × 9 distances (15/25/30/50/100/200/500/1000/2000m)
- `v2_pilgrim_ring_walk_standing.csv` — 36 positions × both real rings
- `v2_elevation_sweep_standing.csv` — vertical directivity, 16 azimuths ×
  13 elevation angles
- `v2_crowd_height_effect.csv` — per-point crowd-adjusted comparison
- `v2_window_driving_conditions.csv` — raw interior driving pressure per
  window per frequency
- `v2_validation.md`, `v2_validation_modes.png` — the analytic mode
  validation gate
- `v2_diff_vs_v1.md` — full round-1-vs-round-2 comparison with all caveats
- `v2_headline_frequency_response.png`, `v2_advantage_zone_directional.png`,
  `v2_height_mechanism_diagram.png` — the corresponding v2 plots

All CSVs have a `#`-prefixed provenance header (geometry version, wall
thickness, window width, producing script) before the data rows — strip
lines starting with `#` before parsing as CSV, same convention throughout.
