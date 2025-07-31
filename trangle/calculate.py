from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, Superimposer, PDBIO
import math
import json
from types import SimpleNamespace
from collections import namedtuple

# === Paths ===
path = Path("/workspaces/Graphormer/TRangle/data")
data_path = path

# Load PCA axes and coreset
with open(data_path / 'coresets.json') as f:
    coresets = json.load(f)

pcL = np.loadtxt(data_path / 'principal_components_alpha.csv')   # Chain A
pcH = np.loadtxt(data_path / 'principal_components_beta.csv')    # Chain B









def save_pymol_visualization_script(
    fname, HL_angle_deg, L_points, H_points, out_script="vis.py"
):
    def cgo_arrow(start, end, radius=0.5, color=(1, 1, 1)):
        return f"""[
    cgo.CYLINDER, {start[0]:.3f}, {start[1]:.3f}, {start[2]:.3f}, {end[0]:.3f}, {end[1]:.3f}, {end[2]:.3f}, {radius},
    {color[0]}, {color[1]}, {color[2]}, {color[0]}, {color[1]}, {color[2]},
    cgo.CONE, {end[0]:.3f}, {end[1]:.3f}, {end[2]:.3f}, {start[0]:.3f}, {start[1]:.3f}, {start[2]:.3f}, {radius*1.6}, 0.0,
    {color[0]}, {color[1]}, {color[2]}, {color[0]}, {color[1]}, {color[2]}, 1.0
]"""

    pdb_name = Path(fname).name
    with open(data_path / out_script, 'w') as f:
        f.write(f"""from pymol import cmd, cgo

cmd.load("{pdb_name}", "structure")
cmd.bg_color('black')
cmd.hide('everything')
cmd.show('cartoon')

L1_vector = {cgo_arrow(L_points.C, L_points.V1, color=(0.0, 0.0, 1.0))}
H1_vector = {cgo_arrow(H_points.C, H_points.V1, color=(0.0, 1.0, 0.0))}
centroid_vector = {cgo_arrow(L_points.C, H_points.C, color=(1.0, 0.5, 0.0))}

cmd.load_cgo(L1_vector, "L1_vector")
cmd.load_cgo(H1_vector, "H1_vector")
cmd.load_cgo(centroid_vector, "centroid_vector")

cmd.pseudoatom(pos={list(L_points.C)}, label="VL_C")
cmd.pseudoatom(pos={list(H_points.C)}, label="VH_C")
cmd.pseudoatom(pos={[((L_points.C[i] + H_points.C[i]) / 2) for i in range(3)]}, label="HL: {HL_angle_deg:.2f} deg")

cmd.zoom("all")
cmd.ray(1200, 1000)
cmd.png("{pdb_name.replace('.pdb', '.png')}", dpi=300)
cmd.save("visualization.pse")
cmd.quit()
""")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python calculate.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    angles, vectors = find_angles(pdb_file, return_vectors=True)
    Lpoints = vectors["Lpoints"]
    Hpoints = vectors["Hpoints"]
    save_pymol_visualization_script("aligned.pdb", angles["HL"], Lpoints, Hpoints)

    print(f"Angles for {pdb_file}:")
    for key, value in angles.items():
        print(f"{key}: {value:.2f} degrees")


