"""
Height effect: does the tower's ~30m elevation help a singer be heard over/
through a gathered crowd, compared to singing from ground level?

REVISION HISTORY (see PROJECT_SUMMARY_FOR_REVIEW.md Sec 4.9 for the full
account -- summarised here since this file has been rebuilt twice in the
same round):

1. Original version borrowed the auditorium "seat-dip effect" (grazing
   attenuation over SEATED audience rows) as an analogy for a standing
   crowd -- always flagged as the weakest empirical input in the project.
2. First revision replaced that with Elena Shabalina's PhD thesis ("The
   Propagation of Low Frequency Sound through an Audience", RWTH Aachen,
   2013), which measures STANDING crowds directly and finds sound decays
   LESS through a crowd than free air at 20-120Hz (impedance-mismatch
   reflection, not absorption). That was a real, correctly-read finding --
   but applying it produced a nonsensical result: a large, confident
   ground-singer ADVANTAGE over the tower, because Shabalina's reflection-
   coefficient formula was extended, ungrounded, across the whole vocal
   range AND applied only to the ground singer's path while the tower's
   path (a geometry her thesis never studied) silently defaulted to zero
   effect. That's not a real finding, it's an artifact of applying an
   evidenced correction to one side of a comparison and leaving the other
   at an unjustified default.
3. THIS version fixes that by going back to basic, frequency-general
   scattering physics instead of extrapolating a narrow low-frequency
   study. A standing person is a scattering obstacle of radius ~0.25m.
   Whether a sound wave "sees" that obstacle is governed by the standard
   Rayleigh-to-geometric scattering transition, controlled by the
   dimensionless parameter ka (k = wavenumber, a = body radius):
     - ka << 1 (wavelength >> body size, roughly Shabalina's own 20-120Hz
       band): bodies are acoustically near-transparent, negligible excess
       loss. This model reproduces that regime as a natural consequence,
       not an assumption -- a genuine cross-check against Shabalina's
       measured result, even though her specific reflection-coefficient
       formula is no longer used.
     - ka >~ 1 (wavelength approaching/below body size -- this starts
       around 200Hz for a 0.25m body, i.e. within this project's own
       practical vocal range): bodies become real acoustic obstacles, and
       a horizontal path through a crowd must cross a sequence of them.
       This is the same physics behind the standard concert-sound
       convention of flying mid/high loudspeakers ("tops") above an
       audience while only bass subwoofers stay grazing at ground level --
       and it only gets more severe climbing toward speech-intelligibility
       frequencies.
   The same physical mechanism (cumulative shadowing along the direction
   of travel) is applied identically to BOTH the ground singer's path
   (long horizontal traverse through the crowd -> large effect) and the
   tower's path (near-vertical incidence, crosses at most a very short
   horizontal stretch of crowd near the listener -> small effect) -- so
   the asymmetry between the two is now a consequence of genuinely
   different geometry under one shared formula, not an artifact of what
   happened to be measured for one side and not the other.

General grounding: this is the standard physics of sound propagating
through a random array of cylindrical scatterers -- the same approach
used for excess attenuation through forests/tree trunks (see e.g. Price,
Attenborough & Heap 1988, cited in Shabalina's own bibliography) and
formally treated by multiple-scattering theory (Twersky 1962/1966, also
in her bibliography). No claim is made here to that literature's precise
coefficients -- this is a simple, clearly-labelled order-of-magnitude
model built on the same well-established qualitative physics, in keeping
with this project's general "broad, honestly-hedged, not falsely precise"
approach.

Still a lightweight geometric + physical post-processing correction
applied on top of the already-computed free-air BEM results (no new BEM
solve of the crowd itself).
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np

BASE = f"{HERE}"

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
H_CROWD = 1.8             # m, assumed standing-pilgrim height -- a guess
MAX_CROWD_DEPTH = 20.0    # m, a bound on how deep a continuously-packed local
                            # gathering plausibly extends between singer and
                            # listener -- NOT a claim about the whole ring
                            # distance (Shabalina's own real crowd depths,
                            # both BEM and live, were in the 10-30m range;
                            # this stops the model from assuming a solid
                            # packed crowd stretches unrealistically far,
                            # e.g. out to the 500-2000m reference distances
                            # used elsewhere in this project)

TOWER_ELEVATION = 30.0
WINDOW_MID_HEIGHT = 0.6
TOWER_SOURCE_HEIGHT = TOWER_ELEVATION + WINDOW_MID_HEIGHT  # ~30.6m
GROUND_SINGER_HEIGHT = 1.75
LISTENER_HEIGHT = 1.6

# ---------------------------------------------------------------------------
# Basic scattering physics: a crowd as a random array of rigid cylinders
# ---------------------------------------------------------------------------
C0 = 340.0           # m/s, speed of sound in air
BODY_RADIUS = 0.25   # m, average human-body radius

NU_PILGRIM = 0.5             # pers/m^2, central "loose gathering" estimate --
                               # anchored to Shabalina's own lowest live-
                               # measured festival-crowd bracket (0.345
                               # pers/m^2), well below her denser BEM/scale
                               # densities (1.3-3.8 pers/m^2, a packed
                               # concert crowd, not a monastic gathering)
NU_SENSITIVITY_RANGE = (0.3, 1.2)  # matches her live-measurement brackets


def scattering_cross_section_per_length(freq):
    """Per-unit-length scattering cross-section of a single infinite rigid
    cylinder (a standing person), as a function of frequency. Standard
    smooth interpolation between the two known asymptotic limits of 2D
    (cylinder) scattering theory: Rayleigh limit sigma_s ~ (ka)^3 at small
    ka, saturating to the geometric-optics limit sigma_s = 2a at large ka
    (the classical "extinction paradox" factor of 2). Not a claim to the
    exact Bessel-function solution -- a simple, standard-shaped, clearly
    labelled order-of-magnitude interpolation."""
    k = 2 * np.pi * freq / C0
    ka = k * BODY_RADIUS
    return 2 * BODY_RADIUS * ka ** 3 / (1 + ka ** 3)


def ka_regime(freq):
    k = 2 * np.pi * freq / C0
    ka = k * BODY_RADIUS
    if ka < 0.3:
        return "rayleigh (bodies near-transparent)"
    if ka < 3.0:
        return "transition"
    return "geometric (bodies are real obstacles)"


def in_layer_path_length(h_source, h_listener, horiz_dist, h_crowd=H_CROWD):
    """Portion of the straight-line source->listener path (in metres of
    horizontal distance) that falls below h_crowd -- the exposure proxy for
    how much of a path crosses a sequence of standing bodies."""
    if h_source <= h_crowd and h_listener <= h_crowd:
        return horiz_dist
    if h_source > h_crowd and h_listener > h_crowd:
        return 0.0
    frac = (h_source - h_crowd) / (h_source - h_listener)
    frac = min(max(frac, 0.0), 1.0)
    return (1 - frac) * horiz_dist


def crowd_effect_db(h_source, h_listener, horiz_dist, freq, nu=NU_PILGRIM):
    """Excess attenuation (dB, always >= 0 -- this model has no "bonus"
    regime) from cumulative scattering by standing bodies along the
    in-layer portion of the source->listener path. SPL_with_crowd =
    SPL_free - crowd_effect_db, so this always reduces level, never
    raises it -- by construction, this model cannot produce a crowd-
    driven advantage for a grazing path over a steep one."""
    L = min(in_layer_path_length(h_source, h_listener, horiz_dist), MAX_CROWD_DEPTH)
    sigma_s = scattering_cross_section_per_length(freq)
    mu = nu * sigma_s  # linear attenuation coefficient, 1/m (Beer-Lambert)
    return 4.343 * mu * L  # 10*log10(e) * mu * L


# ===========================================================================
# Sanity checks: confirm the model reproduces near-zero effect at
# Shabalina's validated low-frequency end (cross-check, not an assumption),
# and confirm the ground singer's penalty grows monotonically and
# substantially through the vocal range while the tower's stays small.
# ===========================================================================
if __name__ == "__main__":
    print("=== Check 1: regime and per-person cross-section vs. frequency ===")
    for f in [50, 80, 120, 165, 200, 250, 300, 390]:
        sigma = scattering_cross_section_per_length(f)
        print(f"  f={f:5.0f}Hz  ka_regime={ka_regime(f):35s}  sigma_s={sigma:.4f} m")
    print("Expect: near-zero sigma_s at 50-80Hz (cross-check vs. Shabalina's "
          "measured near-transparency there); rising toward ~0.5m (2*a) at "
          "high frequency.\n")

    print("=== Check 2: ground singer vs. tower crowd effect across frequency"
          " (d=75m, representative pilgrim-ring distance) ===")
    for f in [50, 80, 120, 165, 196, 250, 300, 390]:
        ge = crowd_effect_db(GROUND_SINGER_HEIGHT, LISTENER_HEIGHT, 75.0, f)
        te = crowd_effect_db(TOWER_SOURCE_HEIGHT, LISTENER_HEIGHT, 75.0, f)
        print(f"  f={f:5.0f}Hz  ground_effect_db={ge:6.2f}  tower_effect_db={te:6.2f}"
              f"  (ground penalty - tower penalty = {ge - te:6.2f})")
    print("Expect: both near zero at low frequency; ground penalty growing "
          "much faster than tower's as frequency rises; tower never worse "
          "off than the ground singer from this mechanism.\n")

# ---------------------------------------------------------------------------
# Apply to the already-computed free-air results
# ---------------------------------------------------------------------------
rows = list(csv.DictReader(open(f"{BASE}/ground_projection_standing.csv")))
for r in rows:
    for k in ("freq_hz", "distance_m", "azimuth_deg", "spl_tower_db", "spl_free_db", "gain_db"):
        r[k] = float(r[k])

# Q used for the free-field baseline is embedded in spl_free_db already
# (same total volume velocity as the tower's actual window output) -- for
# a fair "ground-level singer, same power" comparison we reuse that same
# baseline SPL, just recomputed at a different source height.

out_rows = []
for r in rows:
    d = r["distance_m"]
    f = r["freq_hz"]

    tower_effect = crowd_effect_db(TOWER_SOURCE_HEIGHT, LISTENER_HEIGHT, d, f)
    ground_effect = crowd_effect_db(GROUND_SINGER_HEIGHT, LISTENER_HEIGHT, d, f)

    spl_tower_with_crowd = r["spl_tower_db"] - tower_effect
    spl_ground_with_crowd = r["spl_free_db"] - ground_effect

    out_rows.append({
        "freq_hz": f, "distance_m": d, "azimuth_deg": r["azimuth_deg"],
        "regime": ka_regime(f),
        "tower_crowd_effect_db": tower_effect,
        "ground_crowd_effect_db": ground_effect,
        "spl_tower_with_crowd_db": spl_tower_with_crowd,
        "spl_ground_singer_with_crowd_db": spl_ground_with_crowd,
        "tower_advantage_db": spl_tower_with_crowd - spl_ground_with_crowd,
    })

out_csv = f"{BASE}/crowd_height_effect.csv"
with open(out_csv, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
    writer.writeheader()
    writer.writerows(out_rows)
print(f"wrote {out_csv} ({len(out_rows)} rows)")

# summary
print(f"\n{'dist':>6} {'mean_tower_effect':>18} {'mean_ground_effect':>19} {'mean_tower_advantage_dB':>24}")
for d in sorted(set(r["distance_m"] for r in out_rows)):
    sel = [r for r in out_rows if r["distance_m"] == d]
    ta = np.mean([r["tower_crowd_effect_db"] for r in sel])
    ga = np.mean([r["ground_crowd_effect_db"] for r in sel])
    adv = np.mean([r["tower_advantage_db"] for r in sel])
    print(f"{d:6.0f} {ta:18.2f} {ga:19.2f} {adv:24.2f}")
