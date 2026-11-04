import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import matplotlib.pyplot as plt

BASE = f"{HERE}"


def load(path):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        for k, v in r.items():
            if k != "source_case":
                r[k] = float(v)
    return rows


ground = {sc: load(f"{BASE}/ground_projection_{sc}.csv") for sc in ("standing", "seated")}
elev = {sc: load(f"{BASE}/elevation_sweep_{sc}.csv") for sc in ("standing", "seated")}
driving = list(csv.DictReader(open(f"{BASE}/window_driving_conditions.csv")))
for r in driving:
    for k in ("freq_hz", "window_index", "azimuth_deg", "pressure_re", "pressure_im", "pressure_mag"):
        r[k] = float(r[k])

fig = plt.figure(figsize=(14, 11))

# --- Panel A: gain vs frequency, standing case, one line per distance
ax = fig.add_subplot(2, 2, 1)
freqs = sorted(set(r["freq_hz"] for r in ground["standing"]))
dists = sorted(set(r["distance_m"] for r in ground["standing"]))
colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(dists)))
for d, color in zip(dists, colors):
    means = [np.mean([r["gain_db"] for r in ground["standing"]
                       if r["distance_m"] == d and r["freq_hz"] == f]) for f in freqs]
    ax.plot(freqs, means, marker="o", color=color, label=f"{int(d)} m")
ax.axhline(0, color="k", linewidth=1, linestyle="--")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Tower gain vs. free-field baseline (dB)")
ax.set_title("MAIN TEST: standing, central, omnidirectional source")
ax.legend(title="ground distance", fontsize=8)
ax.grid(alpha=0.3)

# --- Panel B: same for seated control case
ax = fig.add_subplot(2, 2, 2)
for d, color in zip(dists, colors):
    means = [np.mean([r["gain_db"] for r in ground["seated"]
                       if r["distance_m"] == d and r["freq_hz"] == f]) for f in freqs]
    ax.plot(freqs, means, marker="o", color=color, label=f"{int(d)} m")
ax.axhline(0, color="k", linewidth=1, linestyle="--")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Tower gain vs. free-field baseline (dB)")
ax.set_title("CONTROL: seated, aimed at window 0\n(cf. bell tests: known to show no effect)")
ax.legend(title="ground distance", fontsize=8)
ax.grid(alpha=0.3)

# --- Panel C: elevation-angle directivity, standing vs seated, avg over azimuth, at 143Hz
ax = fig.add_subplot(2, 2, 3)
target_f = 143.0
for sc, color, style in [("standing", "orangered", "-"), ("seated", "steelblue", "--")]:
    els = sorted(set(r["elevation_deg"] for r in elev[sc]))
    means = [np.mean([r["gain_db"] for r in elev[sc]
                       if r["freq_hz"] == target_f and r["elevation_deg"] == e]) for e in els]
    ax.plot(els, means, marker="o", color=color, linestyle=style, label=sc)
ax.axvline(0, color="gray", linewidth=1, linestyle=":")
ax.axhline(0, color="k", linewidth=1, linestyle="--")
ax.set_xlabel("Elevation angle (deg): -90=straight down, 0=horizontal, +90=straight up")
ax.set_ylabel("Tower gain vs. free-field baseline (dB)")
ax.set_title(f"Vertical directivity at {int(target_f)} Hz, {200} m range\n(avg over compass azimuth)")
ax.legend()
ax.grid(alpha=0.3)

# --- Panel D: per-window driving pressure, standing vs seated, at 143Hz
ax = fig.add_subplot(2, 2, 4, projection="polar")
for sc, color in [("standing", "orangered"), ("seated", "steelblue")]:
    sel = [r for r in driving if r["source_case"] == sc and r["freq_hz"] == target_f]
    sel.sort(key=lambda r: r["window_index"])
    az = np.radians([r["azimuth_deg"] for r in sel] + [sel[0]["azimuth_deg"] + 360])
    mag = [r["pressure_mag"] for r in sel] + [sel[0]["pressure_mag"]]
    mag_db = 20 * np.log10(np.array(mag) / 20e-6 / np.sqrt(2))
    ax.plot(az, mag_db, marker="o", color=color, label=sc)
ax.set_theta_zero_location("E")
ax.set_theta_direction(1)
ax.set_title(f"Per-window interior driving SPL at {int(target_f)} Hz\n"
             f"(input to the exterior solve, before radiation)")
ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1))

plt.tight_layout()
out = f"{BASE}/radiation_results_summary.png"
plt.savefig(out, dpi=140)
print("saved", out)
