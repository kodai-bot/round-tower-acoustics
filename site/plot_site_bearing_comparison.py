"""
Illustrative (NOT precision) comparison of tower-to-church bearings: this
project's own modeled convention (due west) vs. the three real, quantified
sites from Herity (1984). Distances and bearings only -- no enclosure
outlines, since the source doesn't give them at that level of detail.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import matplotlib.pyplot as plt

OUT = f"{HERE}/site_bearing_comparison.png"

# (site, distance_m, bearing_deg_from_church_to_tower [0=E,90=N,180=W,270=S
#  in standard math convention where 0deg=+x/East, ccw positive], is_corrected)
SITES = [
    ("This project's model\n(due west convention)", 17.0, 180, False),
    ("Clonmacnois", 30.0, 135, False),   # NW of church
    ("Kells", 25.0, 225, True),          # SW of church (corrected from "S")
    ("Glendalough", 15.0, 135, False),   # NW of church
]

fig, axes = plt.subplots(1, 4, figsize=(16, 4.5), subplot_kw={"projection": "polar"})

for ax, (name, dist, bearing_deg, corrected) in zip(axes, SITES):
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)

    # church at centre
    ax.plot(0, 0, marker="+", color="black", markersize=14, markeredgewidth=2)
    ax.annotate("church", xy=(0, 0), xytext=(0.35, 0.15), fontsize=8, ha="center")

    # tower position
    theta = np.radians(bearing_deg)
    ax.plot(theta, dist, marker="s", color="firebrick" if corrected else "seagreen", markersize=12, zorder=5)

    # tower's 4 cardinal windows, drawn as short arrows from the tower position
    for wangle in [0, 90, 180, 270]:
        wx = dist * np.cos(theta) + 5 * np.cos(np.radians(wangle))
        wy = dist * np.sin(theta) + 5 * np.sin(np.radians(wangle))
        w_theta = np.arctan2(wy, wx)
        w_r = np.hypot(wx, wy)
        ax.plot([theta, w_theta], [dist, w_r], color="gray", linewidth=1, alpha=0.6)

    # line from church to tower (the "doorway targets the church" line)
    ax.plot([0, theta], [0, dist], color="black", linewidth=1, linestyle=":", alpha=0.7)

    label = f"{name}\n{dist:.0f}m, bearing {bearing_deg}°"
    if corrected:
        label += "\n(bearing corrected from\nHerity's 'S', see text)"
    ax.set_title(label, fontsize=9)
    ax.set_ylim(0, 35)
    ax.set_rticks([10, 20, 30])
    ax.tick_params(labelsize=6)

fig.suptitle("Tower-to-church bearing: this project's convention vs. 3 real sites (Herity 1984)\n"
             "Illustrative only -- bearings/distances from the source text, NOT digitized survey plans. "
             "Gray lines = the tower's 4 cardinal windows (assumed, not site-confirmed).",
             fontsize=10)
plt.tight_layout(rect=[0, 0, 1, 0.88])
plt.savefig(OUT, dpi=140)
print("saved", OUT)
