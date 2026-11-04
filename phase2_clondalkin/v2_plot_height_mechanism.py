"""v2 height-mechanism explainer diagram -- thin re-pointing of v1's
plot_height_mechanism.py, updated only for v2's tower source height
(window-band midpoint sits higher: 1.765m vs v1's 0.6m -- see
v2_crowd_height_effect.py). Purely illustrative, no CSV data dependency,
same as v1's version."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import v2_params as P
from v2_crowd_height_effect import TOWER_SOURCE_HEIGHT, GROUND_SINGER_HEIGHT, LISTENER_HEIGHT, H_CROWD

OUT = f"{P.MODEL_DIR}/v2_height_mechanism_diagram.png"
LISTENER_DIST = 60.0

fig = plt.figure(figsize=(13, 8.5))
gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1], hspace=0.35)
ax_full = fig.add_subplot(gs[0])
ax_zoom = fig.add_subplot(gs[1])

TOWER_X = 15.0
GROUND_X = -12.0

ax_full.axhline(0, color="saddlebrown", linewidth=2, zorder=1)
ax_full.add_patch(Rectangle((-18, 0), LISTENER_DIST + TOWER_X + 33, H_CROWD, facecolor="firebrick", alpha=0.15, zorder=0))
ax_full.annotate("crowd layer (~1.8m)", xy=(LISTENER_DIST + TOWER_X - 5, H_CROWD), xytext=(LISTENER_DIST + TOWER_X - 5, 4.5),
                  fontsize=9, color="firebrick", ha="center")

ax_full.add_patch(Rectangle((TOWER_X - 1.8, 0), 3.6, TOWER_SOURCE_HEIGHT, facecolor="lightgray", edgecolor="dimgray", zorder=2))
ax_full.plot(TOWER_X, TOWER_SOURCE_HEIGHT, marker="*", color="seagreen", markersize=16, zorder=6)
ax_full.annotate(f"tower singer (~{TOWER_SOURCE_HEIGHT:.1f}m)", xy=(TOWER_X, TOWER_SOURCE_HEIGHT),
                  xytext=(TOWER_X + 7, TOWER_SOURCE_HEIGHT - 1.5),
                  fontsize=9.5, color="seagreen", ha="left", weight="bold")

ax_full.plot(GROUND_X, GROUND_SINGER_HEIGHT, marker="*", color="firebrick", markersize=16, zorder=6)
ax_full.annotate("ground-level singer (1.75m)", xy=(GROUND_X, GROUND_SINGER_HEIGHT), xytext=(GROUND_X, -4.5),
                  fontsize=9.5, color="firebrick", ha="center", weight="bold")

listener_x = TOWER_X + LISTENER_DIST
ax_full.plot(listener_x, LISTENER_HEIGHT, marker="o", color="steelblue", markersize=8, zorder=6)
ax_full.annotate("listener (1.6m)", xy=(listener_x, LISTENER_HEIGHT), xytext=(listener_x + 3, 7),
                  fontsize=9.5, color="steelblue")

ax_full.plot([TOWER_X, listener_x], [TOWER_SOURCE_HEIGHT, LISTENER_HEIGHT], color="seagreen", linewidth=2.2, zorder=5,
             label="tower singer's path: steep, barely touches the crowd layer")
ax_full.plot([GROUND_X, listener_x], [GROUND_SINGER_HEIGHT, LISTENER_HEIGHT], color="firebrick", linewidth=2.2, zorder=5,
             label="ground singer's path: embedded in the crowd the whole way")

ax_full.set_xlim(-22, listener_x + 20)
ax_full.set_ylim(-7, TOWER_SOURCE_HEIGHT + 6)
ax_full.set_ylabel("height (m)")
ax_full.set_title("v2 (Clondalkin dimensions): basic scattering physics -- bodies are obstacles once "
                   "wavelength approaches body size.\nThe ground singer's path crosses a long sequence "
                   "of standing bodies; the tower's steep path\ncrosses almost none -- so the ground "
                   "singer's penalty grows with frequency, the tower's barely does.")
ax_full.legend(loc="upper left", fontsize=9)
ax_full.set_xticks([])

zoom_x0, zoom_x1 = listener_x - 16, listener_x + 4
ax_full.add_patch(Rectangle((zoom_x0, -0.5), zoom_x1 - zoom_x0, 6, fill=False, edgecolor="black",
                             linestyle=":", linewidth=1.2))

ax_zoom.axhline(0, color="saddlebrown", linewidth=2, zorder=1)
ax_zoom.add_patch(Rectangle((zoom_x0, 0), zoom_x1 - zoom_x0, H_CROWD, facecolor="firebrick", alpha=0.15, zorder=0))

xs = np.linspace(zoom_x0, zoom_x1, 100)
tower_ray_y = TOWER_SOURCE_HEIGHT + (LISTENER_HEIGHT - TOWER_SOURCE_HEIGHT) * ((xs - TOWER_X) / (listener_x - TOWER_X))
ground_ray_y = GROUND_SINGER_HEIGHT + (LISTENER_HEIGHT - GROUND_SINGER_HEIGHT) * ((xs - GROUND_X) / (listener_x - GROUND_X))

ax_zoom.plot(xs, tower_ray_y, color="seagreen", linewidth=2.5,
             label="tower singer's path -- only the last stretch dips below crowd height\n(crosses at most one body's worth of obstacle)")
ax_zoom.plot(xs, ground_ray_y, color="firebrick", linewidth=2.5,
             label="ground singer's path -- inside the crowd layer for its entire length\n(crosses a full sequence of bodies -- full penalty applies)")
ax_zoom.plot(listener_x, LISTENER_HEIGHT, marker="o", color="steelblue", markersize=8, zorder=6)

ax_zoom.set_xlim(zoom_x0, zoom_x1)
ax_zoom.set_ylim(-0.5, 5.5)
ax_zoom.set_xlabel("horizontal distance -- zoomed detail near the listener")
ax_zoom.set_ylabel("height (m)")
ax_zoom.legend(loc="upper left", fontsize=7.5)
ax_zoom.set_title("Zoomed: this is why the penalty applies almost entirely to the ground singer's path", fontsize=10)

plt.tight_layout()
plt.savefig(OUT, dpi=140)
print("saved", OUT)
