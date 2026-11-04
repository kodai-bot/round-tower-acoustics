import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import matplotlib.pyplot as plt

BASE = f"{HERE}"

rows = list(csv.DictReader(open(f"{BASE}/pilgrim_ring_walk_standing.csv")))
for r in rows:
    for k in ("freq_hz", "theta_site_deg", "dist_from_tower_m", "spl_tower_db", "spl_free_db", "gain_db"):
        r[k] = float(r[k])

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), subplot_kw={"projection": "polar"})

show_freqs = [50, 100, 143, 175, 250]
colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(show_freqs)))

for ax, ring in zip(axes, ("2nd_ring", "3rd_ring")):
    for f, color in zip(show_freqs, colors):
        sel = sorted([r for r in rows if r["ring"] == ring and r["freq_hz"] == f],
                     key=lambda r: r["theta_site_deg"])
        theta = np.radians([r["theta_site_deg"] for r in sel] + [sel[0]["theta_site_deg"] + 360])
        spl = [r["spl_tower_db"] for r in sel] + [sel[0]["spl_tower_db"]]
        ax.plot(theta, spl, color=color, label=f"{int(f)} Hz", linewidth=1.6)
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    # mark tower side (west, theta=180) and church-facing side (east, theta=0)
    ax.annotate("nearest tower\n(west side)", xy=(np.pi, ax.get_rmax()), ha="center", fontsize=8, color="dimgray")
    ax.annotate("church-facing side\n(farthest from tower)", xy=(0, ax.get_rmax()), ha="center", fontsize=8, color="dimgray")
    ax.set_title(f"{ring.replace('_',' ')} -- SPL (dB) while circumambulating")

axes[0].legend(loc="upper left", bbox_to_anchor=(-0.35, 1.15), fontsize=8)

plt.tight_layout()
out = f"{BASE}/pilgrim_ring_walk.png"
plt.savefig(out, dpi=140)
print("saved", out)
