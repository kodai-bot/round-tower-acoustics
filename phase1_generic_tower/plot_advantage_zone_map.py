import os
HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import sys
sys.path.insert(0, f"{HERE}")
from crowd_height_effect import crowd_effect_db, TOWER_SOURCE_HEIGHT, GROUND_SINGER_HEIGHT, LISTENER_HEIGHT
import csv
import numpy as np

BASE = f"{HERE}"
OUT = f"{BASE}/advantage_zone_map.png"

TOWER_OFFSET = 17.0
RING_1, RING_2, RING_3 = 20.0, 45.0, 91.5
TARGET_FREQ = 165.0  # representative frequency within the sean-nos tonic
                      # range -- NOT the old seat-dip peak (150Hz), which no
                      # longer has any special status; see crowd_height_effect.py
                      # for the basic scattering-physics model now used

# gain_db (structural, from the actual BEM solve) vs distance -- interpolate
# from the existing ground-projection results at the representative frequency
rows = list(csv.DictReader(open(f"{BASE}/ground_projection_standing.csv")))
for r in rows:
    for k in ("freq_hz", "distance_m", "gain_db"):
        r[k] = float(r[k])
freqs = sorted(set(r["freq_hz"] for r in rows))
nearest_f = min(freqs, key=lambda x: abs(x - TARGET_FREQ))
dist_gain = sorted({(r["distance_m"], np.mean([rr["gain_db"] for rr in rows
                     if rr["freq_hz"] == nearest_f and rr["distance_m"] == r["distance_m"]]))
                     for r in rows if r["freq_hz"] == nearest_f})
d_grid = np.array([d for d, g in dist_gain])
g_grid = np.array([g for d, g in dist_gain])


def structural_gain_db(d):
    return np.interp(d, d_grid, g_grid)


def net_advantage_db(d):
    te = crowd_effect_db(TOWER_SOURCE_HEIGHT, LISTENER_HEIGHT, d, nearest_f)
    ge = crowd_effect_db(GROUND_SINGER_HEIGHT, LISTENER_HEIGHT, d, nearest_f)
    return structural_gain_db(d) - te + ge  # tower(with crowd effect) - ground(with
    # crowd effect); relative gain already embeds the free-field comparison for both.
    # NOTE: crowd_effect_db is now always >= 0 (a basic scattering/shadowing PENALTY,
    # never a bonus -- see crowd_height_effect.py) and grows with frequency, applying
    # much more strongly to the ground singer's long horizontal path through the crowd
    # than to the tower's near-vertical one. This typically makes net_advantage_db
    # MORE positive (tower more favoured) than the structural gain alone, growing with
    # frequency -- the opposite of the intermediate Shabalina-only revision.


fig, (ax_map, ax_curve) = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={"width_ratios": [1.15, 1]})

# ---- left: site plan colored by net advantage as a function of tower distance ----
radii = np.linspace(1, 130, 300)
advantages = np.array([net_advantage_db(r) for r in radii])

theta = np.linspace(0, 2 * np.pi, 200)
R, TH = np.meshgrid(radii, theta)
ADV = np.tile(advantages, (len(theta), 1))
X = TOWER_OFFSET * -1 + R * np.cos(TH)  # centered on TOWER position, not the site origin,
Y = R * np.sin(TH)                       # since advantage depends on distance FROM TOWER

vmax = np.abs(advantages).max()
mesh = ax_map.pcolormesh(X, Y, ADV, cmap="RdYlGn", vmin=-vmax, vmax=vmax, shading="auto")
cb = fig.colorbar(mesh, ax=ax_map, label=f"tower advantage over ground singer at {nearest_f:.0f}Hz (dB)")

# overlay the monastery rings (centred on the SITE origin, offset from the tower)
for r, label in [(RING_1, "1st ring"), (RING_2, "2nd ring"), (RING_3, "3rd ring")]:
    ax_map.add_patch(Circle((0, 0), r, fill=False, edgecolor="black", linewidth=1.3))
    ax_map.annotate(label, xy=(0, r), xytext=(2, r + 2), fontsize=8.5)

ax_map.plot(-TOWER_OFFSET, 0, marker="s", color="black", markersize=8)
ax_map.annotate("tower", xy=(-TOWER_OFFSET, 0), xytext=(-TOWER_OFFSET - 4, -8), fontsize=9, weight="bold")
ax_map.plot(0, 0, marker="+", color="black", markersize=10)

ax_map.set_xlim(-135, 105)
ax_map.set_ylim(-120, 120)
ax_map.set_aspect("equal")
ax_map.set_xlabel("metres")
ax_map.set_ylabel("metres")
ax_map.set_title(f"Net tower advantage by distance (basic scattering-physics crowd model)\n(green = tower wins, red = ground singer wins, at {nearest_f:.0f}Hz)")

# ---- right: the advantage curve itself, with ring distances marked ----
ax_curve.plot(radii, advantages, color="black", linewidth=2)
ax_curve.axhline(0, color="gray", linewidth=1, linestyle="--")
ax_curve.fill_between(radii, advantages, 0, where=advantages > 0, color="seagreen", alpha=0.25)
ax_curve.fill_between(radii, advantages, 0, where=advantages <= 0, color="firebrick", alpha=0.25)

ymin, ymax = min(advantages.min(), -1) - 1, max(advantages.max(), 1) + 1
ax_curve.set_ylim(ymin, ymax)
for r, label, color in [(28, "2nd ring (near)", "steelblue"), (62, "2nd ring (far)", "steelblue"),
                         (74.5, "3rd ring (near)", "darkorange"), (108.5, "3rd ring (far)", "darkorange")]:
    ax_curve.axvline(r, color=color, linewidth=1, linestyle=":")
    ax_curve.annotate(label, xy=(r, ymin + 0.5), xytext=(r, ymin + 0.5),
                       fontsize=7.5, color=color, ha="center", rotation=90, va="bottom")

ax_curve.set_xlabel("distance from tower (m)")
ax_curve.set_ylabel("net tower advantage over ground singer (dB)")
ax_curve.set_title(f"Advantage vs. distance at {nearest_f:.0f}Hz\n(basic scattering physics: tower favoured, growing with frequency -- see PROJECT_SUMMARY_FOR_REVIEW.md Sec 4.9)")
ax_curve.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUT, dpi=140)
print("saved", OUT)
