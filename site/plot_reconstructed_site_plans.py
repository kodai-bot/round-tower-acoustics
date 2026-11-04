"""
Reconstructed (approximate, prose-derived) schematic plans for the three
Herity (1984) sites with round towers and quantified monument positions:
Clonmacnois, Kells, Glendalough. NOT digitized survey plans -- positions
are computed from the paper's stated distances/bearings from the church,
each with its exact source quote. Where only a direction is given (no
distance), the monument is drawn as a ray from the church, not a fixed
point, and clearly marked as such.

Cross-to-tower distances are DERIVED (triangulated from two church-relative
vectors), not directly stated by Herity -- flagged in the output.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

OUT = f"{HERE}/reconstructed_site_plans.png"


def compass_xy(dist, compass_deg):
    math_deg = 90 - compass_deg
    r = np.radians(math_deg)
    return dist * np.cos(r), dist * np.sin(r)


# (name, distance_or_None, compass_bearing_deg, marker, color, source_quote)
SITES = {
    "Clonmacnois": {
        "scale": 45,
        "points": [
            ("church (Daimhliag Mor)", 0, None, "+", "black", "reference point"),
            ("South Cross", np.hypot(4, 12), np.degrees(np.arctan2(4, -12)) % 360, "P", "darkorange",
             "'about 12m south and 4m east of its present facade'"),
            ("North Cross", np.hypot(4, 12), np.degrees(np.arctan2(4, 12)) % 360, "P", "darkorange",
             "'the same distance north, and similarly placed east'"),
            ("Round Tower", 30, 315, "s", "seagreen", "'30m north-west of the west doorway'"),
            ("Cross of the Scriptures", 12.5, 270, "P", "firebrick",
             "'12.5m from the west facade' (direction: 'this open space', assumed west)"),
        ],
        "rays": [
            ("Teampall Ciaran", 45, "gray", "'north-east of the Daimhliag Mor' (no distance given)"),
        ],
    },
    "Kells": {
        "scale": 100,
        "points": [
            ("church (present church)", 0, None, "+", "black", "reference point"),
            ("Round Tower", 25, 225, "s", "seagreen",
             "'to the south, 25m away' -- CORRECTED to SW, see text"),
            ("Cross of Patrick & Columba", 25, 225, "P", "firebrick",
             "'beside it' (the tower) -- co-located, not independently measured"),
            ("St Columb's House", 90, 315, "D", "steelblue", "'north-west of the west facade, 90m away'"),
        ],
        "rays": [
            ("Broken Cross", 270, "gray", "'in the open space opposite the west doorway' (no distance given)"),
        ],
    },
    "Glendalough": {
        "scale": 40,
        "points": [
            ("church (cathedral)", 0, None, "+", "black", "reference point"),
            ("Round Tower", 15, 315, "s", "seagreen", "'15m away' to the NW of the facade"),
            ("Priest's House", 27, 202.5, "D", "steelblue",
             "'27m to the south-south-west' (on a lower terrace)"),
        ],
        "rays": [
            ("plain granite High Cross", 180, "gray", "'south of the west facade' (no distance given)"),
        ],
    },
}

fig, axes = plt.subplots(1, 3, figsize=(17, 6.5))

# manual label offsets (as fraction of panel scale) for points that sit too
# close together to label cleanly at a fixed small offset
LABEL_OFFSETS = {
    "Cross of the Scriptures": (0.03, -0.10),
    "church (Daimhliag Mor)": (0.03, 0.08),
    "Cross of Patrick & Columba": (0.05, -0.14),
    "Round Tower": (0.03, 0.05),
}

for ax, (site, data) in zip(axes, SITES.items()):
    scale = data["scale"]
    tower_xy = None
    for name, dist, bearing, marker, color, note in data["points"]:
        if bearing is None:
            x, y = 0, 0
        else:
            x, y = compass_xy(dist, bearing)
        if "Round Tower" in name:
            tower_xy = (x, y)
        ax.scatter([x], [y], marker=marker, color=color, s=120, zorder=5, edgecolor="black", linewidth=0.6)
        dx, dy = LABEL_OFFSETS.get(name, (0.03, 0.03))
        ax.annotate(name, xy=(x, y), xytext=(x + scale * dx, y + scale * dy), fontsize=7.5)

    # direction-only rays
    for name, bearing, color, note in data["rays"]:
        x, y = compass_xy(scale * 0.55, bearing)
        ax.annotate("", xy=(x, y), xytext=(0, 0),
                     arrowprops=dict(arrowstyle="->", color=color, linewidth=1.3, linestyle="--"))
        ax.annotate(f"{name}\n(direction only)", xy=(x, y), xytext=(x, y), fontsize=7, color=color, ha="center")

    ax.plot(0, 0, marker="+", color="black", markersize=0)  # keep origin in view
    ax.set_xlim(-scale, scale)
    ax.set_ylim(-scale, scale)
    ax.set_aspect("equal")
    ax.axhline(0, color="lightgray", linewidth=0.5, zorder=0)
    ax.axvline(0, color="lightgray", linewidth=0.5, zorder=0)
    ax.set_title(site, fontsize=12, weight="bold")
    ax.set_xlabel("metres (E +)")
    if site == "Clonmacnois":
        ax.set_ylabel("metres (N +)")

fig.suptitle("Reconstructed (approximate) core layouts from Herity (1984) -- prose-derived, NOT digitized survey plans\n"
             "Solid markers = distance+bearing both stated; dashed rays = bearing only, distance not given",
             fontsize=10)
plt.tight_layout(rect=[0, 0, 1, 0.90])
plt.savefig(OUT, dpi=140)
print("saved", OUT)

# print the derived tower-monument distances table
print("\n=== Derived (triangulated, NOT directly stated) distances from tower ===")
for site, data in SITES.items():
    pts = {}
    for name, dist, bearing, marker, color, note in data["points"]:
        pts[name] = (0, 0) if bearing is None else compass_xy(dist, bearing)
    tower = next((v for k, v in pts.items() if "Round Tower" in k), None)
    print(f"\n{site}:")
    for name, (x, y) in pts.items():
        if tower is None or name == next(k for k in pts if "Round Tower" in k):
            continue
        d = np.hypot(x - tower[0], y - tower[1])
        print(f"  {name:28s} -> tower: {d:6.2f}m (derived)")
