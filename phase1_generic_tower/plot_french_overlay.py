import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = f"{HERE}"

rows = list(csv.DictReader(open(f"{BASE}/ground_projection_standing.csv")))
for r in rows:
    for k in ("freq_hz", "distance_m", "azimuth_deg", "spl_tower_db", "gain_db"):
        r[k] = float(r[k])

freqs = sorted(set(r["freq_hz"] for r in rows))
gain_100m = [np.mean([r["gain_db"] for r in rows if r["freq_hz"] == f and r["distance_m"] == 100]) for f in freqs]

fig, ax = plt.subplots(figsize=(13, 6.5))
ax.plot(freqs, gain_100m, color="darkorange", linewidth=1.8, zorder=2, label="model: tower gain vs free-field @100m")
ax.axhline(0, color="k", linewidth=1, linestyle="--")

# French (major key) scale-degree time profile, "not scientific" per the
# source -- Fifth+Fourth given as a combined ~25% "launchpad" figure and
# Second+Sixth as a combined <15% "passing" figure, split evenly here
# since no finer breakdown was given.
notes = {
    "Do tonic": (0, 35, "steelblue"),
    "Re passing": (2, 7.5, "gray"),
    "Mi major 3rd": (4, 20, "goldenrod"),
    "Fa 4th (launchpad)": (5, 12.5, "seagreen"),
    "Sol 5th (launchpad,\n'camped out' on)": (7, 12.5, "crimson"),
    "La passing": (9, 7.5, "gray"),
}

for tonic_hz, marker, tonic_label in [(196.0, "o", "G3 tonic"), (207.65, "s", "Ab3 tonic")]:
    for note, (semis, pct, color) in notes.items():
        freq = tonic_hz * 2 ** (semis / 12)
        size = 60 + pct * 14
        g = np.interp(freq, freqs, gain_100m)
        ax.scatter([freq], [g], s=size, marker=marker, color=color,
                   edgecolor="black", linewidth=0.6, zorder=5, alpha=0.9)

note_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markeredgecolor="black",
                        markersize=10, label=f"{n.splitlines()[0]} ({pct}%)")
                for n, (s, pct, c) in notes.items()]
tonic_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor="dimgray", markeredgecolor="black",
                         markersize=8, label="if tonic=G3 (196Hz)"),
                 Line2D([0], [0], marker="s", color="w", markerfacecolor="dimgray", markeredgecolor="black",
                        markersize=8, label="if tonic=Ab3 (207.65Hz)")]
leg1 = ax.legend(handles=note_handles, loc="lower left", title="French scale degree\n(marker size ~ % of time spent)", fontsize=8)
ax.add_artist(leg1)
ax.legend(handles=tonic_handles, loc="upper right", fontsize=8)

ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Tower gain vs. free-field baseline (dB) @ 100m")
ax.set_title("Where the notes of 'French' fall on the tower's projection curve\n"
              "(mesh now validated to 390Hz -- Sol/5th and La/6th are real computed points, not extrapolated)")
ax.grid(alpha=0.3)
ax.set_xlim(140, 360)

plt.tight_layout()
out = f"{BASE}/french_overlay.png"
plt.savefig(out, dpi=140)
print("saved", out)
