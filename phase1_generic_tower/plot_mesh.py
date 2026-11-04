import os
HERE = os.path.dirname(os.path.abspath(__file__))

import gmsh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

gmsh.initialize()
gmsh.open(f"{HERE}/tower_upper_chamber.msh")

node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
coord_map = {t: node_coords[3*i:3*i+3] for i, t in enumerate(node_tags)}

fig = plt.figure(figsize=(10, 9))
ax = fig.add_subplot(111, projection='3d')

colors = {"windows": "orangered", "walls": "lightsteelblue"}

for dim, tag in gmsh.model.getEntities(2):
    phys = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
    if len(phys) == 0:
        continue
    name = gmsh.model.getPhysicalName(dim, phys[0])
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(dim, tag)
    for et, ents, ent_nodes in zip(elem_types, elem_tags, elem_node_tags):
        props = gmsh.model.mesh.getElementProperties(et)
        nnodes = props[3]
        ent_nodes = np.array(ent_nodes).reshape(-1, nnodes)
        tris = [[coord_map[n] for n in tri] for tri in ent_nodes]
        pc = Poly3DCollection(tris, facecolor=colors.get(name, "gray"),
                               edgecolor='k', linewidths=0.1, alpha=0.9)
        ax.add_collection3d(pc)

ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_zlim(0, 3.6)
ax.set_box_aspect([1,1,1])
ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
ax.set_title("Round tower upper chamber -- exterior BEM skin\n(orange = window/Dirichlet, blue = wall+dome/Neumann)")
ax.view_init(elev=22, azim=35)

plt.tight_layout()
out = f"{HERE}/tower_mesh_preview.png"
plt.savefig(out, dpi=140)
print("saved", out)
gmsh.finalize()
