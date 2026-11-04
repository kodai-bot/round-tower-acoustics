import os
HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Arc

OUT = f"{HERE}/tower_cross_section.png"

TOWER_RADIUS = 2.2
WALL_HEIGHT = 2.5
DOME_RISE = 0.9
WINDOW_HEIGHT = 1.45
WINDOW_WIDTH = 0.46
TOWER_ELEVATION = 30.0
SOURCE_STANDING_Z = 1.75
SOURCE_SEATED_Z = 0.7
LISTENER_HEIGHT = 1.6

fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={"width_ratios": [1, 1.6]})

# --- left: full elevation view, tower on its plinth above the ground
ax_full.add_patch(Rectangle((-TOWER_RADIUS, 0), 2 * TOWER_RADIUS, TOWER_ELEVATION, facecolor="lightgray", edgecolor="dimgray"))
ax_full.add_patch(Rectangle((-TOWER_RADIUS, TOWER_ELEVATION), 2 * TOWER_RADIUS, WALL_HEIGHT, facecolor="wheat", edgecolor="black"))
theta = np.linspace(0, np.pi, 100)
dome_x = TOWER_RADIUS * np.cos(theta)
dome_y = TOWER_ELEVATION + WALL_HEIGHT + DOME_RISE * np.sin(theta)
ax_full.plot(dome_x, dome_y, color="black")
ax_full.fill_between(dome_x, TOWER_ELEVATION + WALL_HEIGHT, dome_y, color="wheat", edgecolor="black")
ax_full.axhline(0, color="saddlebrown", linewidth=2)
ax_full.annotate("ground level", xy=(6, 0), fontsize=8, color="saddlebrown")
ax_full.annotate(f"~{TOWER_ELEVATION:.0f}m\n(90-110ft estimate,\nmidpoint used)", xy=(-6.5, TOWER_ELEVATION / 2),
                  fontsize=8, ha="center", color="dimgray")
ax_full.plot([-5.5, -5.5], [0, TOWER_ELEVATION], color="gray", linewidth=0.8)
ax_full.plot([-6, -5], [0, 0], color="gray", linewidth=0.8)
ax_full.plot([-6, -5], [TOWER_ELEVATION, TOWER_ELEVATION], color="gray", linewidth=0.8)

listener_x = 30
ax_full.plot(listener_x, LISTENER_HEIGHT, marker="o", color="steelblue", markersize=6)
ax_full.annotate("pilgrim\n(1.6m)", xy=(listener_x, LISTENER_HEIGHT), xytext=(listener_x + 2, 6), fontsize=8, color="steelblue")

ax_full.set_xlim(-12, 45)
ax_full.set_ylim(-3, 40)
ax_full.set_aspect("equal")
ax_full.set_title("Tower elevation above ground")
ax_full.set_xlabel("metres")
ax_full.set_ylabel("metres")
ax_full.grid(alpha=0.2)

# --- right: zoomed upper chamber with source positions and window band
ax_zoom.add_patch(Rectangle((-TOWER_RADIUS, 0), 2 * TOWER_RADIUS, WALL_HEIGHT, facecolor="wheat", edgecolor="black", alpha=0.4))
dome_y2 = WALL_HEIGHT + DOME_RISE * np.sin(theta)
ax_zoom.plot(dome_x, dome_y2, color="black")
ax_zoom.fill_between(dome_x, WALL_HEIGHT, dome_y2, color="wheat", edgecolor="black", alpha=0.4)

# window band (front + back, since we're viewing a vertical slice through window 0/window 2 axis)
for side in (-1, 1):
    ax_zoom.add_patch(Rectangle((side * TOWER_RADIUS - (0.15 if side > 0 else 0), 0), 0.15, WINDOW_HEIGHT,
                                 facecolor="orangered", edgecolor="black", alpha=0.85, zorder=5))
ax_zoom.annotate(f"window\n(0-{WINDOW_HEIGHT:.2f}m,\n{WINDOW_WIDTH:.2f}m wide)", xy=(TOWER_RADIUS, WINDOW_HEIGHT / 2),
                  xytext=(TOWER_RADIUS + 0.6, WINDOW_HEIGHT / 2), fontsize=8, color="orangered")

ax_zoom.plot(0, SOURCE_STANDING_Z, marker="*", color="darkorange", markersize=18, zorder=6)
ax_zoom.annotate(f"standing source\n({SOURCE_STANDING_Z}m) -- main test", xy=(0, SOURCE_STANDING_Z),
                  xytext=(0.3, SOURCE_STANDING_Z + 0.5), fontsize=8, color="darkorange")

ax_zoom.plot(TOWER_RADIUS - 0.3, SOURCE_SEATED_Z, marker="*", color="steelblue", markersize=18, zorder=6)
ax_zoom.annotate(f"seated source\n({SOURCE_SEATED_Z}m) -- control", xy=(TOWER_RADIUS - 0.3, SOURCE_SEATED_Z),
                  xytext=(-2.1, SOURCE_SEATED_Z - 0.6), fontsize=8, color="steelblue")

ax_zoom.axhline(WINDOW_HEIGHT, color="gray", linewidth=0.6, linestyle=":")
ax_zoom.annotate("standing mouth height (1.75m) is\nwell ABOVE the window tops (1.45m) --\nno direct line-of-sight aiming",
                  xy=(-2.1, 1.9), fontsize=8, color="dimgray")

ax_zoom.set_xlim(-3.2, 3.5)
ax_zoom.set_ylim(-0.3, 4.0)
ax_zoom.set_aspect("equal")
ax_zoom.set_title("Upper chamber detail\n(radius 2.2m, wall height 2.5m, dome rise 0.9m)")
ax_zoom.set_xlabel("metres")
ax_zoom.grid(alpha=0.2)

plt.tight_layout()
plt.savefig(OUT, dpi=140)
print("saved", OUT)
