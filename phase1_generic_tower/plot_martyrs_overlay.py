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
gain_100m = [np.mean([r["gain_db"] for r in rows if r["freq_hz"] == f and r["distance_m"] == 100]) for f in freqs]

fig, ax = plt.subplots(figsize=(12, 6.5))
ax.plot(freqs, gain_100m, color="darkorange", linewidth=1.8, zorder=2, label="model: tower gain vs free-field @100m")
ax.axhline(0, color="k", linewidth=1, linestyle="--")

# Martyrs scale-degree time profile, transposed onto the precentor's
# actual sung tonic (not the tune's nominal "D") -- two scenarios for
# where the tonic actually sits (G3 vs Ab3)
notes = {
    "D tonic": (0, 37, "steelblue"),
    "E passing": (2, 8, "gray"),
    "F minor 3rd\n(emotional core)": (3, 25, "crimson"),
    "G 4th": (5, 17, "seagreen"),
    "C subtonic": (-2, 8, "gray"),
}

for tonic_hz, marker, tonic_label in [(196.0, "o", "G3 tonic"), (207.65, "s", "Ab3 tonic")]:
    for note, (semis, pct, color) in notes.items():
        freq = tonic_hz * 2 ** (semis / 12)
        # interpolate model gain at this exact frequency
        g = np.interp(freq, freqs, gain_100m)
        size = 60 + pct * 14
        ax.scatter([freq], [g], s=size, marker=marker, color=color,
                   edgecolor="black", linewidth=0.6, zorder=5, alpha=0.85)

# legend proxies
from matplotlib.lines import Line2D
note_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markeredgecolor="black",
                        markersize=10, label=f"{n.splitlines()[0]} ({pct}%)")
                for n, (s, pct, c) in notes.items()]
tonic_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor="dimgray", markeredgecolor="black",
                         markersize=8, label="if tonic=G3 (196Hz)"),
                 Line2D([0], [0], marker="s", color="w", markerfacecolor="dimgray", markeredgecolor="black",
                        markersize=8, label="if tonic=Ab3 (207.65Hz)")]
leg1 = ax.legend(handles=note_handles, loc="lower left", title="Martyrs scale degree\n(marker size ~ % of time spent)", fontsize=8)
ax.add_artist(leg1)
ax.legend(handles=tonic_handles, loc="lower right", fontsize=8)

ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Tower gain vs. free-field baseline (dB) @ 100m")
ax.set_title("Where the notes of 'Martyrs' actually fall on the tower's projection curve\n"
              "(D=tonic/'safe harbour' sits near the model's best gain; F=minor 3rd/'emotional core', held and ornamented, sits near its worst)")
ax.grid(alpha=0.3)
ax.set_xlim(140, 300)

plt.tight_layout()
out = f"{BASE}/martyrs_overlay.png"
plt.savefig(out, dpi=140)
print("saved", out)
