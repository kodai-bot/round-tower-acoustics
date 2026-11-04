import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import matplotlib.pyplot as plt

BASE = f"{HERE}"

rows = list(csv.DictReader(open(f"{BASE}/ground_projection_standing.csv")))
for r in rows:
    for k in ("freq_hz", "distance_m", "azimuth_deg", "spl_tower_db", "gain_db"):
        r[k] = float(r[k])

freqs = sorted(set(r["freq_hz"] for r in rows))
dists = [50, 100, 200, 500]

fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

ax = axes[0]
colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(dists)))
for d, color in zip(dists, colors):
    spl = [np.mean([r["spl_tower_db"] for r in rows if r["freq_hz"] == f and r["distance_m"] == d]) for f in freqs]
    ax.plot(freqs, spl, marker="o", markersize=3, color=color, label=f"{int(d)}m")
ax.set_ylabel("SPL (dB)")
ax.set_title("Fine-resolution frequency response (standing/central source) vs. Gaelic psalm precentor range")
ax.legend(title="distance")
ax.grid(alpha=0.3)

ax2 = axes[1]
gain = [np.mean([r["gain_db"] for r in rows if r["freq_hz"] == f and r["distance_m"] == 100]) for f in freqs]
ax2.plot(freqs, gain, marker="o", markersize=3, color="darkorange")
ax2.axhline(0, color="k", linewidth=1, linestyle="--")
ax2.set_ylabel("Gain vs. free-field (dB) @ 100m")
ax2.set_xlabel("Frequency (Hz)")
ax2.grid(alpha=0.3)

# mark precentor range on both panels
for ax in axes:
    ax.axvspan(165, 311, color="tab:red", alpha=0.06)
    ax.axvspan(196, 207, color="tab:red", alpha=0.15)
    ax.axvline(196, color="firebrick", linestyle=":", linewidth=1)
    ax.axvline(207, color="firebrick", linestyle=":", linewidth=1)

axes[0].annotate("precentor full range\n(Eb3-Eb4, ~165-311Hz)", xy=(238, axes[0].get_ylim()[1]*0.95),
                  fontsize=8, color="firebrick", ha="center", va="top")
axes[0].annotate("tonic\n(G3-Ab3)", xy=(201, axes[0].get_ylim()[1]*0.7),
                  fontsize=8, color="firebrick", ha="center")

plt.tight_layout()
out = f"{BASE}/psalm_range_comparison.png"
plt.savefig(out, dpi=140)
print("saved", out)
