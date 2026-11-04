"""
Refined version of plot_advantage_zone_map.py: uses the actual azimuth-
resolved gain_db (not azimuthally averaged) to check whether the "net
tower advantage" zone has real directional structure -- e.g. lobed toward
or away from the window directions -- rather than assuming it's a uniform
circle around the tower. The crowd-clearing correction itself doesn't
depend on azimuth (it's a pure distance/height effect), so any directional
skew in the result comes entirely from the tower's own interference
pattern (real, confirmed structure -- see the 250Hz on-axis/diagonal
findings).
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import sys
sys.path.insert(0, f"{HERE}")
from crowd_height_effect import crowd_effect_db, TOWER_SOURCE_HEIGHT, GROUND_SINGER_HEIGHT, LISTENER_HEIGHT

BASE = f"{HERE}"
OUT = f"{BASE}/advantage_zone_directional.png"

TOWER_OFFSET = 17.0
RING_1, RING_2, RING_3 = 20.0, 45.0, 91.5
TARGET_FREQ = 165.0  # sean-nos tonic range, see crowd_height_effect.py

rows = list(csv.DictReader(open(f"{BASE}/ground_projection_standing.csv")))
for r in rows:
    for k in ("freq_hz", "distance_m", "azimuth_deg", "gain_db"):
        r[k] = float(r[k])
freqs = sorted(set(r["freq_hz"] for r in rows))
nearest_f = min(freqs, key=lambda x: abs(x - TARGET_FREQ))

dists = sorted(set(r["distance_m"] for r in rows))
azs = sorted(set(r["azimuth_deg"] for r in rows))
gain_lookup = {(r["distance_m"], r["azimuth_deg"]): r["gain_db"]
               for r in rows if r["freq_hz"] == nearest_f}


def structural_gain_db(d, az_deg):
    az_deg = az_deg % 360
    nearest_az = min(azs, key=lambda a: abs((a - az_deg + 180) % 360 - 180))
    d_arr = np.array(dists)
    g_arr = np.array([gain_lookup[(dd, nearest_az)] for dd in dists])
    return np.interp(d, d_arr, g_arr)


def net_advantage_db(d, az_deg):
    te = crowd_effect_db(TOWER_SOURCE_HEIGHT, LISTENER_HEIGHT, d, nearest_f)
    ge = crowd_effect_db(GROUND_SINGER_HEIGHT, LISTENER_HEIGHT, d, nearest_f)
    return structural_gain_db(d, az_deg) - te + ge


fig, ax_map = plt.subplots(figsize=(8.5, 8))

radii = np.linspace(1, 130, 220)
theta_deg = np.linspace(0, 360, 180)
R, TH = np.meshgrid(radii, theta_deg)
ADV = np.vectorize(net_advantage_db)(R, TH)

X = -TOWER_OFFSET + R * np.cos(np.radians(TH))
Y = R * np.sin(np.radians(TH))

vmax = np.abs(ADV).max()
mesh = ax_map.pcolormesh(X, Y, ADV, cmap="RdYlGn", vmin=-vmax, vmax=vmax, shading="auto")
fig.colorbar(mesh, ax=ax_map, label=f"tower advantage over ground singer at {nearest_f:.0f}Hz (dB)")

for r, label in [(RING_1, "1st ring"), (RING_2, "2nd ring"), (RING_3, "3rd ring")]:
    ax_map.add_patch(Circle((0, 0), r, fill=False, edgecolor="black", linewidth=1.3))
    ax_map.annotate(label, xy=(0, r), xytext=(2, r + 2), fontsize=8.5)

# mark the 4 window directions (window 0 faces +x / east, toward the church)
for i, ang in enumerate([0, 90, 180, 270]):
    dx = 25 * np.cos(np.radians(ang))
    dy = 25 * np.sin(np.radians(ang))
    ax_map.annotate("", xy=(-TOWER_OFFSET + dx, dy), xytext=(-TOWER_OFFSET, 0),
                     arrowprops=dict(arrowstyle="->", color="black", linewidth=1.2))
ax_map.annotate("window 0\n(faces church)", xy=(-TOWER_OFFSET + 25, 0), xytext=(-TOWER_OFFSET + 28, -8), fontsize=7.5)

ax_map.plot(-TOWER_OFFSET, 0, marker="s", color="black", markersize=8)
ax_map.annotate("tower", xy=(-TOWER_OFFSET, 0), xytext=(-TOWER_OFFSET - 4, -12), fontsize=9, weight="bold")
ax_map.plot(0, 0, marker="+", color="black", markersize=10)

ax_map.set_xlim(-135, 105)
ax_map.set_ylim(-120, 120)
ax_map.set_aspect("equal")
ax_map.set_xlabel("metres")
ax_map.set_ylabel("metres")
ax_map.set_title(f"Directionally-resolved advantage zone at {nearest_f:.0f}Hz (basic scattering-physics crowd model)\n"
                  f"(NOT azimuth-averaged -- shows the tower's real 4-fold interference structure)")

plt.tight_layout()
plt.savefig(OUT, dpi=140)
print("saved", OUT)
