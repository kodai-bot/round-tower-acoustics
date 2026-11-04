import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, f"{HERE}")
from crowd_height_effect import crowd_effect_db, TOWER_SOURCE_HEIGHT, GROUND_SINGER_HEIGHT, LISTENER_HEIGHT

BASE = f"{HERE}"

rows = list(csv.DictReader(open(f"{BASE}/ground_projection_standing.csv")))
for r in rows:
    for k in ("freq_hz", "distance_m", "spl_tower_db", "spl_free_db"):
        r[k] = float(r[k])

freqs = sorted(set(r["freq_hz"] for r in rows))
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

ax = axes[0]
for d, color in zip([50, 100, 200, 500], plt.cm.viridis(np.linspace(0.1, 0.85, 4))):
    adv = []
    for f in freqs:
        sel = [r for r in rows if r["freq_hz"] == f and r["distance_m"] == d]
        te = crowd_effect_db(TOWER_SOURCE_HEIGHT, LISTENER_HEIGHT, d, f)
        ge = crowd_effect_db(GROUND_SINGER_HEIGHT, LISTENER_HEIGHT, d, f)
        spl_t = np.mean([r["spl_tower_db"] for r in sel]) - te
        spl_g = np.mean([r["spl_free_db"] for r in sel]) - ge
        adv.append(spl_t - spl_g)
    ax.plot(freqs, adv, marker="o", markersize=3, color=color, label=f"{d}m")
ax.axhline(0, color="k", linewidth=1, linestyle="--")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Tower advantage over ground-level singer (dB)\n(structural diffraction + basic-physics scattering/shadowing crowd model)")
ax.set_title("Net height effect vs. frequency (rebuilt on basic scattering physics)\n(positive = tower wins; crowd penalty grows with frequency as bodies\nstop being acoustically transparent -- see PROJECT_SUMMARY_FOR_REVIEW.md Sec 4.9)")
ax.legend(title="distance")
ax.grid(alpha=0.3)

ax2 = axes[1]
ring_rows = list(csv.DictReader(open(f"{BASE}/pilgrim_ring_walk_standing.csv")))
for r in ring_rows:
    for k in ("freq_hz", "dist_from_tower_m", "spl_tower_db", "spl_free_db"):
        r[k] = float(r[k])

for ring, color in [("2nd_ring", "firebrick"), ("3rd_ring", "steelblue")]:
    rfreqs = sorted(set(r["freq_hz"] for r in ring_rows))
    adv = []
    for f in rfreqs:
        sel = [r for r in ring_rows if r["ring"] == ring and r["freq_hz"] == f]
        vals = []
        for r in sel:
            d = r["dist_from_tower_m"]
            te = crowd_effect_db(TOWER_SOURCE_HEIGHT, LISTENER_HEIGHT, d, f)
            ge = crowd_effect_db(GROUND_SINGER_HEIGHT, LISTENER_HEIGHT, d, f)
            vals.append((r["spl_tower_db"] - te) - (r["spl_free_db"] - ge))
        adv.append(np.mean(vals))
    ax2.plot(rfreqs, adv, marker="o", markersize=3, color=color, label=ring.replace("_", " "))
ax2.axhline(0, color="k", linewidth=1, linestyle="--")
ax2.axvspan(131, 196, color="purple", alpha=0.08, label="sean-nos tonic range (131-196Hz)")
ax2.set_xlabel("Frequency (Hz)")
ax2.set_ylabel("Tower advantage over ground-level singer (dB)")
ax2.set_title("At actual pilgrim ring positions (rebuilt model)\n(28-62m / 74.5-108.5m from tower)")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)

plt.tight_layout()
out = f"{BASE}/crowd_height_effect.png"
plt.savefig(out, dpi=140)
print("saved", out)
