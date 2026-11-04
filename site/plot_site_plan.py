import os
HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

OUT = f"{HERE}/site_plan.png"

TOWER_OFFSET = 17.0
RING_1 = 20.0   # inner sanctum / clergy
RING_2 = 45.0   # pilgrims, circumambulation
RING_3 = 91.5   # outer boundary
TOWER_RADIUS = 2.2

fig, ax = plt.subplots(figsize=(9, 9))

for r, label, style in [(RING_1, "1st ring (clergy / church)", "-"),
                         (RING_2, "2nd ring (pilgrims)", "-"),
                         (RING_3, "3rd ring (outer boundary)", "-")]:
    ax.add_patch(Circle((0, 0), r, fill=False, edgecolor="steelblue", linewidth=1.4, linestyle=style))
    ax.annotate(label, xy=(0, r), xytext=(3, r + 2), fontsize=9, color="steelblue")

ax.plot(0, 0, marker="+", color="black", markersize=12, markeredgewidth=2)
ax.annotate("common origin\n(church / tabernacle)", xy=(0, 0), xytext=(4, -6), fontsize=8, color="dimgray")

tower_x, tower_y = -TOWER_OFFSET, 0.0
ax.add_patch(Circle((tower_x, tower_y), TOWER_RADIUS, facecolor="firebrick", edgecolor="black", alpha=0.85, zorder=5))
ax.annotate("tower", xy=(tower_x, tower_y), xytext=(tower_x - 8, tower_y + 6), fontsize=10, color="firebrick", weight="bold")

# 4 cardinal windows, window 0 facing +x (east, toward the church)
for i, ang in enumerate([0, 90, 180, 270]):
    dx = 4.5 * np.cos(np.radians(ang))
    dy = 4.5 * np.sin(np.radians(ang))
    ax.annotate("", xy=(tower_x + dx, tower_y + dy), xytext=(tower_x, tower_y),
                arrowprops=dict(arrowstyle="->", color="orangered", linewidth=1.5))
ax.annotate("window 0 faces\nthe church doorway", xy=(tower_x + 4.5, tower_y), xytext=(tower_x + 8, tower_y - 8),
            fontsize=8, color="orangered", arrowprops=dict(arrowstyle="-", color="orangered", linewidth=0.6))

# mark near/far side distances on 2nd ring
near_x = -RING_2
far_x = RING_2
ax.plot([tower_x, near_x], [0, 0], color="gray", linewidth=0.8, linestyle=":")
ax.plot([tower_x, far_x], [0, 0], color="gray", linewidth=0.8, linestyle=":")
ax.annotate(f"{abs(near_x - tower_x):.0f}m", xy=((tower_x + near_x) / 2, 2), fontsize=8, color="gray", ha="center")
ax.annotate(f"{abs(far_x - tower_x):.0f}m", xy=((tower_x + far_x) / 2, 2), fontsize=8, color="gray", ha="center")

ax.set_xlim(-100, 100)
ax.set_ylim(-100, 100)
ax.set_aspect("equal")
ax.set_xlabel("metres (east +)")
ax.set_ylabel("metres (north +)")
ax.set_title("Monastic enclosure site plan\n(tower sits at the inner edge of the 2nd ring, 17m west of common origin)")
ax.grid(alpha=0.2)

import matplotlib.lines as mlines
plt.tight_layout()
plt.savefig(OUT, dpi=140)
print("saved", OUT)
