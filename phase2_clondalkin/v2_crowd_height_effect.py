"""
v2 crowd-clearing / height-effect post-processing -- re-pointing of v1's
crowd_height_effect.py (its "Attempt 3" cylinder-scattering physics,
described in full there) at v2's new BEM results, per
rebuild_spec.md SS6 step 7 ("height comparison").

No logic changes from v1's model -- the §8 carry-forward item this
defers is specifically the FINITE-IMPEDANCE-BOUNDARY upgrade (Shabalina
equivalent-fluid parameters as a real BEM boundary condition), not this
existing lightweight post-processing correction, which stays in scope as
part of the core SS6 pipeline. See v1's crowd_height_effect.py docstring
for the full physics rationale (Rayleigh-to-geometric cylinder-scattering
transition, applied identically to both the tower's near-vertical path
and a ground singer's long horizontal path through the crowd) -- not
repeated here in full, only the v2-specific parameter change is:

  - WINDOW_MID_HEIGHT: v1 used 0.6m (~half of its 1.2192m window height,
    with the window band starting at floor level). v2's window band sits
    HIGHER and is TALLER (sill 1.219m, head 2.311m, per the RIA drawing --
    see v2_params.py), so the source-height proxy used for "where does the
    sound plausibly exit the tower" is v2_params.WINDOW_SILL_Z to
    WINDOW_HEAD_Z's midpoint, not v1's value.
"""
import csv
import numpy as np

import v2_params as P

# ---------------------------------------------------------------------
H_CROWD = 1.8
MAX_CROWD_DEPTH = 20.0

WINDOW_MID_HEIGHT = (P.WINDOW_SILL_Z + P.WINDOW_HEAD_Z) / 2.0   # 1.765m
TOWER_SOURCE_HEIGHT = P.TOWER_ELEVATION + WINDOW_MID_HEIGHT
GROUND_SINGER_HEIGHT = 1.75
LISTENER_HEIGHT = P.LISTENER_HEIGHT

C0 = 340.0
BODY_RADIUS = 0.25
NU_PILGRIM = 0.5


def scattering_cross_section_per_length(freq):
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
    if h_source <= h_crowd and h_listener <= h_crowd:
        return horiz_dist
    if h_source > h_crowd and h_listener > h_crowd:
        return 0.0
    frac = (h_source - h_crowd) / (h_source - h_listener)
    frac = min(max(frac, 0.0), 1.0)
    return (1 - frac) * horiz_dist


def crowd_effect_db(h_source, h_listener, horiz_dist, freq, nu=NU_PILGRIM):
    L = min(in_layer_path_length(h_source, h_listener, horiz_dist), MAX_CROWD_DEPTH)
    sigma_s = scattering_cross_section_per_length(freq)
    mu = nu * sigma_s
    return 4.343 * mu * L


if __name__ == "__main__":
    print(f"v2 tower source height: {TOWER_SOURCE_HEIGHT:.2f}m "
          f"(elevation {P.TOWER_ELEVATION}m + window-band midpoint {WINDOW_MID_HEIGHT:.3f}m)")

    with open(f"{P.MODEL_DIR}/v2_ground_projection_standing.csv") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    rows = list(csv.DictReader(lines))
    for r in rows:
        for k in ("freq_hz", "distance_m", "azimuth_deg", "spl_tower_db", "spl_free_db", "gain_db"):
            r[k] = float(r[k])

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

    out_csv = f"{P.MODEL_DIR}/v2_crowd_height_effect.csv"
    with open(out_csv, "w", newline="") as fh:
        for line in P.provenance_header("v2_crowd_height_effect.py"):
            fh.write(line + "\n")
        writer = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"wrote {out_csv} ({len(out_rows)} rows)")

    print(f"\n{'dist':>6} {'mean_tower_effect':>18} {'mean_ground_effect':>19} {'mean_tower_advantage_dB':>24}")
    for d in sorted(set(r["distance_m"] for r in out_rows)):
        sel = [r for r in out_rows if r["distance_m"] == d]
        ta = np.mean([r["tower_crowd_effect_db"] for r in sel])
        ga = np.mean([r["ground_crowd_effect_db"] for r in sel])
        adv = np.mean([r["tower_advantage_db"] for r in sel])
        print(f"{d:6.0f} {ta:18.2f} {ga:19.2f} {adv:24.2f}")
