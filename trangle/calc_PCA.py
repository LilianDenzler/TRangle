from sklearn.decomposition import PCA
import os
import pathlib
import json
import numpy as np
from Bio.PDB import PDBParser


def compute_pca_frame(coords):
    """Returns principal component axes from coreset CA coordinates."""
    pca = PCA(n_components=3)
    pca.fit(coords)
    origin = coords.mean(axis=0)
    axes = pca.components_  # rows are PC1, PC2, PC3import numpy as np
from Bio.PDB import PDBParser

def get_ca_coords(pdb_file, chain_id, coreset_residues):
    """Extract CA atom coordinates of coreset residues from PDB."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("consensus", pdb_file)
    model = structure[0]
    chain = model[chain_id]
    coords = []
    for res_id in coreset_residues:
        try:
            res = chain[(" ", res_id, " ")]
            coords.append(res["CA"].get_coord())
        except KeyError:
            print(f"Missing CA for residue {res_id}")
            continue
    return np.array(coords)

def compute_pcs(coords):
    """Compute principal components (eigenvectors of covariance matrix)."""
    coords_centered = coords - coords.mean(axis=0)
    u, s, vh = np.linalg.svd(coords_centered)
    return vh  # shape (3, 3); rows = PC1, PC2, PC3

def save_pca_frame_pymol_script(pdb_path, pcs, origin, out_script="pca_frame_vis.py", label="VL"):
    from pathlib import Path
    import numpy as np

    pdb_name = Path(pdb_path).name

    scale = 10.0  # Adjust arrow length for visibility

    axes = {
        'X': pcs[0],
        'Y': pcs[1],
        'Z': pcs[2]
    }
    colors = {
        'X': (1.0, 0.0, 0.0),  # red
        'Y': (0.0, 1.0, 0.0),  # green
        'Z': (0.0, 0.0, 1.0),  # blue
    }

    def cgo_arrow(start, end, color):
        return f"""[
    cgo.CYLINDER, {start[0]}, {start[1]}, {start[2]}, {end[0]}, {end[1]}, {end[2]}, 0.4,
    {color[0]}, {color[1]}, {color[2]}, {color[0]}, {color[1]}, {color[2]},
    cgo.CONE, {end[0]}, {end[1]}, {end[2]}, {start[0]}, {start[1]}, {start[2]}, 0.8, 0.0,
    {color[0]}, {color[1]}, {color[2]}, {color[0]}, {color[1]}, {color[2]}, 1.0
]"""

    with open(out_script, "w") as f:
        f.write(f"""from pymol import cmd, cgo

cmd.load("{pdb_path}", "structure")
cmd.bg_color("black")
cmd.hide("everything")
cmd.show("cartoon")

""")

        for axis, vec in axes.items():
            start = origin
            end = origin + vec * scale
            f.write(f"{axis}_arrow = {cgo_arrow(start, end, colors[axis])}\n")
            f.write(f'cmd.load_cgo({axis}_arrow, "{label}_{axis}_axis")\n')

        f.write(f"""
cmd.pseudoatom(pos={list(origin)}, label="{label}_origin")
cmd.zoom("all")
cmd.ray(1200, 1000)
cmd.png("{label}_pca_frame.png", dpi=300)
cmd.quit()
""")


# === INPUT ===
path = pathlib.Path(__file__).parent
data_path = path.parent/'data'

# load in coreset residue number dictionary
with open(data_path/'coresets.json') as f:
    coresets = json.load(f)
core_alpha = coresets['A']
core_beta  = coresets['B']

alpha_pdb = data_path/"consensus_A.pdb"
beta_pdb  = data_path/"consensus_B.pdb"

# === RUN ===
coords_A = get_ca_coords(alpha_pdb, "A", core_alpha)
coords_B = get_ca_coords(beta_pdb, "B", core_beta)

pcs_A = compute_pcs(coords_A)  # → save as principal_components_light.csv
pcs_B = compute_pcs(coords_B)  # → save as principal_components_heavy.csv
originA=coords_A.mean(axis=0)
originB=coords_B.mean(axis=0)

save_pca_frame_pymol_script(alpha_pdb, pcs_A, originA, out_script="pca_alpha_vis.py", label="VL")
save_pca_frame_pymol_script(beta_pdb, pcs_B, originB, out_script="pca_beta_vis.py", label="VH")

# Save them
np.savetxt(data_path/"principal_components_alpha.csv", pcs_A)
np.savetxt(data_path/"principal_components_beta.csv", pcs_B)

