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

fig, ax = plt.subplots(figsize=(12, 5.5))
ax.plot(freqs, gain_100m, color="black", linewidth=2)
ax.axhline(0, color="gray", linewidth=1, linestyle="--")
ax.fill_between(freqs, gain_100m, -16, where=[g > -3 for g in gain_100m], color="seagreen", alpha=0.15)

ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Tower gain vs. free-field baseline (dB) @ 100m")
ax.set_title("Tower exterior projection efficiency vs. frequency\n"
              "-- computed purely from room geometry, before any vocal/singing data was introduced --")
ax.grid(alpha=0.3)
ax.set_xlim(45, 395)
ax.set_ylim(-16, 2)

plt.tight_layout()
out = f"{BASE}/headline_frequency_response.png"
plt.savefig(out, dpi=140)
print("saved", out)
