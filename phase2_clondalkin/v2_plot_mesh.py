"""Render the two phase-2 BEM meshes side by side: the exterior skin (what the
outdoor air sees) and the interior chamber (Variant A, ceiled).
Orange = window openings, blue = stone wall/roof. Writes v2_mesh_preview.png."""
import os

import gmsh
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
COLOURS = {"windows": "orangered", "walls": "lightsteelblue"}


def draw(ax, path, title, cut_front=False, wall_alpha=0.95):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(path)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = dict(zip(tags, coords.reshape(-1, 3)))
    for dim, tag in gmsh.model.getEntities(2):
        phys = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
        if not len(phys):
            continue
        name = gmsh.model.getPhysicalName(dim, phys[0])
        for et, _, nodes in zip(*gmsh.model.mesh.getElements(dim, tag)):
            n = gmsh.model.mesh.getElementProperties(et)[3]
            tris = np.array([[xyz[k] for k in tri] for tri in np.asarray(nodes).reshape(-1, n)])
            if cut_front:  # remove the near half so the inside is visible
                tris = tris[tris[:, :, 0].mean(axis=1) + tris[:, :, 1].mean(axis=1) < 0.2]
            ax.add_collection3d(Poly3DCollection(tris, facecolor=COLOURS.get(name, "gray"),
                                                 edgecolor="k", linewidths=0.08,
                                                 alpha=0.95 if name == "windows" else wall_alpha))
    gmsh.finalize()
    ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.set_zlim(0, 5)
    ax.set_box_aspect([4, 4, 5])
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
    ax.view_init(elev=20, azim=35)
    ax.set_title(title)


fig = plt.figure(figsize=(13, 7))
draw(fig.add_subplot(121, projection="3d"), f"{HERE}/v2_mesh_exterior.msh",
     "Exterior skin (wall shown translucent): 3.47 m across,\ncone roof, windows as ducts through a 0.745 m wall", wall_alpha=0.25)
draw(fig.add_subplot(122, projection="3d"), f"{HERE}/v2_mesh_interior_A.msh",
     "Interior chamber (cut away): 1.98 m across, 2.90 m high,\nceiled at the corbel line", cut_front=True)
fig.suptitle("Phase 2 BEM meshes, Clondalkin bell storey (orange = window openings, blue = stone)", y=0.98)
fig.subplots_adjust(top=0.86, left=0.02, right=0.98, wspace=0.05)
out = f"{HERE}/v2_mesh_preview.png"
plt.savefig(out, dpi=140)
print("saved", out)
