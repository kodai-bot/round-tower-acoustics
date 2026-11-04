"""
NOTE: v2_diff_vs_v1.md was edited by hand after generation (the caveat sections).
Re-running this script overwrites those edits; regenerate into a copy and merge.

One-page diff: v1 (superseded) vs v2 (Clondalkin-specific) headline
numbers, per rebuild_spec.md SS7 ("a one-page diff
summarising how headline numbers moved vs. the old (void) results").

Every number here is recomputed directly from each round's own CSVs
(not copied from PROJECT_SUMMARY_FOR_REVIEW.md's prose) so this script
stays correct if either round's data is ever regenerated.

v1's numbers are READ-ONLY inputs -- this script does not modify any v1
file, consistent with the rest of this rebuild.
"""
import csv

import numpy as np

import v2_params as P

V1_GROUND_CSV = f"{P.MODEL_DIR}/../phase1_generic_tower/ground_projection_standing.csv"
V1_DRIVING_CSV = f"{P.MODEL_DIR}/../phase1_generic_tower/window_driving_conditions.csv"
V1_ELEV_CSV = f"{P.MODEL_DIR}/../phase1_generic_tower/elevation_sweep_standing.csv"

V2_GROUND_CSV = f"{P.MODEL_DIR}/v2_ground_projection_standing.csv"
V2_DRIVING_CSV = P.DRIVING_CSV_PATH
V2_ELEV_CSV = f"{P.MODEL_DIR}/v2_elevation_sweep_standing.csv"

OUT_MD = f"{P.MODEL_DIR}/v2_diff_vs_v1.md"


def _read_csv(path, float_cols, skip_hash=False):
    with open(path) as fh:
        lines = [ln for ln in fh if not (skip_hash and ln.startswith("#"))]
    rows = list(csv.DictReader(lines))
    for r in rows:
        for c in float_cols:
            if c in r:
                r[c] = float(r[c])
    return rows


def headline_peak_trough(ground_rows, eval_distance=100.0, exclude=()):
    freqs = sorted(set(r["freq_hz"] for r in ground_rows) - set(exclude))
    gain_by_f = {f: np.mean([r["gain_db"] for r in ground_rows
                              if r["freq_hz"] == f and r["distance_m"] == eval_distance])
                 for f in freqs}
    best_f = max(gain_by_f, key=gain_by_f.get)
    worst_f = min(gain_by_f, key=gain_by_f.get)
    return (best_f, gain_by_f[best_f]), (worst_f, gain_by_f[worst_f])


def window_imbalance(driving_rows):
    freqs = sorted(set(r["freq_hz"] for r in driving_rows))
    imbalances = []
    for f in freqs:
        mags = [r["pressure_mag"] for r in driving_rows if r["freq_hz"] == f]
        imbalances.append(20 * np.log10(max(mags) / max(min(mags), 1e-12)))
    return float(np.mean(imbalances)), float(np.max(imbalances))


def hemisphere_split(elev_rows):
    upper = lower = 0.0
    for row in elev_rows:
        w = np.cos(np.radians(row["elevation_deg"]))
        I = 10 ** (row["spl_tower_db"] / 10)
        if row["elevation_deg"] >= 0:
            upper += I * w
        else:
            lower += I * w
    return lower / (upper + lower)


if __name__ == "__main__":
    v1_ground = _read_csv(V1_GROUND_CSV, ["freq_hz", "distance_m", "azimuth_deg", "gain_db"])
    v1_driving = _read_csv(V1_DRIVING_CSV, ["freq_hz", "pressure_mag"])
    v1_elev = _read_csv(V1_ELEV_CSV, ["elevation_deg", "spl_tower_db"])

    v2_ground = _read_csv(V2_GROUND_CSV, ["freq_hz", "distance_m", "azimuth_deg", "gain_db"], skip_hash=True)
    v2_driving = _read_csv(V2_DRIVING_CSV, ["freq_hz", "pressure_mag"], skip_hash=True)
    v2_elev = _read_csv(V2_ELEV_CSV, ["elevation_deg", "spl_tower_db"], skip_hash=True)

    v1_best, v1_worst = headline_peak_trough(v1_ground)
    v2_best, v2_worst = headline_peak_trough(v2_ground, exclude=P.NONCONVERGED_HZ)  # converged only

    v1_imb_mean, v1_imb_max = window_imbalance(v1_driving)
    v2_imb_mean, v2_imb_max = window_imbalance(v2_driving)

    v1_lower = hemisphere_split(v1_elev)
    v2_lower = hemisphere_split(v2_elev)

    with open(OUT_MD, "w") as fh:
        fh.write("# v2 vs v1 headline diff\n\n")
        fh.write("**All v1 frequency-domain numbers below are SUPERSEDED.** v1 modelled "
                 "the chamber as a zero-thickness shell at a single compromise radius "
                 "(2.2m, a general cross-tower average) -- roughly 2x the true Clondalkin "
                 "internal diameter, with wall thickness omitted entirely and windows as "
                 "flush patches rather than real ducts. v2 rebuilds on RIA-drawing-derived, "
                 "Clondalkin-specific dimensions with a real ~0.745m wall and rectangular "
                 "duct windows. Nothing below should be read as a correction to v1's "
                 "numbers -- it is a different, more physically-grounded model of a "
                 "specific real tower.\n\n")

        fh.write("## Headline projection efficiency (gain vs. free-field @ 100m)\n\n")
        fh.write("| | v1 (superseded, zero-thickness shell) | v2 (Clondalkin, real duct wall) |\n")
        fh.write("|---|---|---|\n")
        fh.write(f"| Best frequency | {v1_best[0]:.0f}Hz, {v1_best[1]:.1f}dB | "
                 f"{v2_best[0]:.0f}Hz, {v2_best[1]:.1f}dB |\n")
        fh.write(f"| Worst frequency | {v1_worst[0]:.0f}Hz, {v1_worst[1]:.1f}dB | "
                 f"{v2_worst[0]:.0f}Hz, {v2_worst[1]:.1f}dB |\n")

        fh.write("\n## Four-window balance (centred standing source)\n\n")
        fh.write("| | v1 | v2 |\n|---|---|---|\n")
        fh.write(f"| Mean imbalance across sweep | {v1_imb_mean:.2f}dB | {v2_imb_mean:.2f}dB |\n")
        fh.write(f"| Worst-case imbalance | {v1_imb_max:.2f}dB | {v2_imb_max:.2f}dB |\n")

        fh.write("\n## Vertical directivity (power-weighted lower-hemisphere fraction)\n\n")
        fh.write("| | v1 | v2 |\n|---|---|---|\n")
        fh.write(f"| Fraction radiated below horizontal | {v1_lower:.1%} | {v2_lower:.1%} |\n")
        fh.write(f"\n(v1 elevation grid: 8 azimuths, unchanged this round in v1 -- see "
                 f"PROJECT_SUMMARY_FOR_REVIEW.md SS5 item 1. v2 elevation grid: "
                 f"{len(P.ELEVATION_AZIMUTHS_DEG)} azimuths.)\n")

        fh.write("\n## What changed in the geometry itself\n\n")
        fh.write("| Parameter | v1 | v2 |\n|---|---|---|\n")
        fh.write("| Interior radius | 2.2m (general average) | 0.9905m (Clondalkin RIA drawing) |\n")
        fh.write(f"| Wall thickness | 0 (zero-thickness shell) | {P.WALL_THICKNESS}m (real, ducted) |\n")
        fh.write("| Window model | flush Dirichlet patch | rectangular duct through full wall depth |\n")
        fh.write(f"| Chamber height | 1.2192m (cylinder) + 2.5908m cone | {P.CHAMBER_HEIGHT}m "
                 f"(Variant A, ceiled -- cone excluded from acoustic domain) |\n")
        fh.write(f"| Source position | (0,0,1.75m), above window band | "
                 f"{P.SOURCE_POSITION}, inside window band |\n")
        fh.write(f"| Frequency ceiling | 390Hz | {max(P.FREQUENCIES):.0f}Hz |\n")
        fh.write(f"| Ground azimuth resolution | 36 points (widened mid-round after an aliasing bug) | "
                 f"{len(P.AZIMUTHS_DEG)} points (from the start) |\n")
        fh.write("| Mesh topology | open tube, missing floor cap (found during v2 planning, "
                 "not corrected in v1 -- see v2_check_watertight.py) | closed, watertight, "
                 "checked automatically |\n")

    print(f"wrote {OUT_MD}")
    with open(OUT_MD) as fh:
        print(fh.read())
