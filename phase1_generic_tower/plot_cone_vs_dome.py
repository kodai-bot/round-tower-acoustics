import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import matplotlib.pyplot as plt

BASE = f"{HERE}"

# current file = cone (just solved); we don't have a separately-saved dome
# CSV anymore since it was overwritten, but we captured the dome-era PNG.
# Reconstruct the dome curve isn't possible without rerunning -- instead
# just plot the current (cone) curve clearly labeled, referencing the
# saved old PNG for comparison. (Kept simple deliberately.)

rows = list(csv.DictReader(open(f"{BASE}/ground_projection_standing.csv")))
for r in rows:
    for k in ("freq_hz", "distance_m", "spl_tower_db", "gain_db"):
        r[k] = float(r[k])

freqs = sorted(set(r["freq_hz"] for r in rows))
gain_100m = [np.mean([r["gain_db"] for r in rows if r["freq_hz"] == f and r["distance_m"] == 100]) for f in freqs]

fig, ax = plt.subplots(figsize=(12, 5.5))
ax.plot(freqs, gain_100m, color="firebrick", linewidth=2, label="CONE roof (Clondalkin scan)")
ax.axhline(0, color="gray", linewidth=1, linestyle="--")
ax.axvspan(130, 196, color="steelblue", alpha=0.12, label="sean-nos tonic range (C3-G3)")

ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Tower gain vs. free-field baseline (dB) @ 100m")
ax.set_title("Cone-roof geometry (Clondalkin scan) -- exterior projection efficiency\n"
              "compare against headline_frequency_response_DOME_OLD.png (shallow-dome version)")
ax.legend()
ax.grid(alpha=0.3)
ax.set_xlim(45, 395)
ax.set_ylim(-16, 2)

plt.tight_layout()
out = f"{BASE}/headline_frequency_response_cone_labeled.png"
plt.savefig(out, dpi=140)
print("saved", out)
