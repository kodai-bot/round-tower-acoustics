"""v2 headline frequency-response plot.

Converged frequencies are drawn as the solid curve; the frequencies in
v2_params.NONCONVERGED_HZ (GMRES hit its iteration cap) are drawn as hollow
grey markers so they can't be read as findings."""
import csv
import numpy as np
import matplotlib.pyplot as plt

import v2_params as P

with open(f"{P.MODEL_DIR}/v2_ground_projection_standing.csv") as fh:
    lines = [ln for ln in fh if not ln.startswith("#")]
rows = list(csv.DictReader(lines))
for r in rows:
    for k in ("freq_hz", "distance_m", "azimuth_deg", "spl_tower_db", "gain_db"):
        r[k] = float(r[k])

freqs = sorted(set(r["freq_hz"] for r in rows))
gain_100m = [np.mean([r["gain_db"] for r in rows if r["freq_hz"] == f and r["distance_m"] == 100]) for f in freqs]

ok = [f not in P.NONCONVERGED_HZ for f in freqs]
f_ok = [f for f, k in zip(freqs, ok) if k]
g_ok = [g for g, k in zip(gain_100m, ok) if k]
f_bad = [f for f, k in zip(freqs, ok) if not k]
g_bad = [g for g, k in zip(gain_100m, ok) if not k]

fig, ax = plt.subplots(figsize=(12, 5.5))
ax.plot(f_ok, g_ok, color="black", linewidth=2, marker="o", markersize=3.5,
        label=f"converged ({len(f_ok)} frequencies)")
ax.plot(f_bad, g_bad, linestyle="none", marker="o", markersize=6, markerfacecolor="none",
        markeredgecolor="gray", label=f"did not converge ({len(f_bad)}) -- not a finding")
ax.axhline(0, color="gray", linewidth=1, linestyle="--")
ax.axvspan(245, max(freqs) + 5, color="gray", alpha=0.07)
ax.legend(loc="upper right")

ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Tower gain vs. free-field baseline (dB) @ 100m")
ax.set_title("v2 (Clondalkin, real wall thickness + duct windows): "
              "exterior projection efficiency vs. frequency\n"
              "-- azimuth-averaged at 100 m, Variant A (ceiled); trustworthy range about 50-330 Hz --")
ax.grid(alpha=0.3)
ax.set_xlim(min(freqs) - 5, max(freqs) + 5)

plt.tight_layout()
out = f"{P.MODEL_DIR}/v2_headline_frequency_response.png"
plt.savefig(out, dpi=140)
print("saved", out)
