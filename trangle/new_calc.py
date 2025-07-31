"""
DESCRIPTION
    Calculate orientation between the alpha and beta domains in TCR Fv regions.

OUTPUT
    - Angles of the Fv regions for each input PDB
    - PyMOL visualisation script and PNG

AUTHOR
    Adapted from ABangle (Dunbar et al., 2012)
"""

import os, json, math
from pathlib import Path
from collections import namedtuple
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Superimposer
from number import write_renumbered_fv

# -------------------------------------------------------
# Data structures
# -------------------------------------------------------
Points = namedtuple("Points", ["C", "V1", "V2"])

# -------------------------------------------------------
# Helpers
# -------------------------------------------------------
def get_coreset_atoms(structure, chain, coresets):
    return [
        atom for atom in structure.get_atoms()
        if atom.get_name() == "CA"
        and atom.parent.parent.id == chain
        and atom.parent.id[1] in coresets[chain]
    ]

def apply_transform(structure, rot, tran):
    for atom in structure.get_atoms():
        atom.transform(rot, tran)

def centroid_and_vectors(structure, chain, pcs, coresets):
    ca_coords = [
        atom.get_coord() for atom in structure.get_atoms()
        if atom.get_name() == "CA"
        and atom.parent.parent.id == chain
        and atom.parent.id[1] in coresets[chain]
    ]
    centroid = np.mean(ca_coords, axis=0)
    coefs, *_ = np.linalg.lstsq(pcs.T, centroid, rcond=None)
    C = coefs @ pcs
    return Points(C=C, V1=C + pcs[0], V2=C + pcs[1])

def as_unit(v): return v / np.linalg.norm(v)
def angle(v1, v2): return math.degrees(math.acos(np.clip(np.dot(v1, v2), -1.0, 1.0)))

def add_cgo_arrow(start,end,color,radius=0.3):
    return f"""[
        cgo.CYLINDER,{start[0]:.3f},{start[1]:.3f},{start[2]:.3f},
                     {end[0]:.3f},{end[1]:.3f},{end[2]:.3f},
                     {radius},
                     {color[0]},{color[1]},{color[2]},
                     {color[0]},{color[1]},{color[2]},
        cgo.CONE,{end[0]:.3f},{end[1]:.3f},{end[2]:.3f},
                 {start[0]:.3f},{start[1]:.3f},{start[2]:.3f},
                 {radius*1.5},0.0,
                 {color[0]},{color[1]},{color[2]},
                 {color[0]},{color[1]},{color[2]},1.0
    ]"""

# -------------------------------------------------------
# Core alignment & angle calculation
# -------------------------------------------------------
def process(input_pdb, consA_pdb, consB_pdb, pcsA, pcsB, coresets, out_dir):
    parser = PDBParser(QUIET=True)
    input_struct = parser.get_structure("input", input_pdb)
    consA = parser.get_structure("consA", consA_pdb)
    consB = parser.get_structure("consB", consB_pdb)

    # Step 1: superpose input A onto consensus A
    si = Superimposer()
    si.set_atoms(get_coreset_atoms(consA,"A",coresets),
                 get_coreset_atoms(input_struct,"A",coresets))
    rotA, tranA = si.rotran
    apply_transform(input_struct, rotA, tranA)

    # Step 2: superpose consensus B onto moved input B
    si2 = Superimposer()
    si2.set_atoms(get_coreset_atoms(input_struct,"B",coresets),
                  get_coreset_atoms(consB,"B",coresets))
    rotB, tranB = si2.rotran
    apply_transform(consB, rotB, tranB)

    Apts = centroid_and_vectors(consA,"A",pcsA,coresets)
    Bpts = centroid_and_vectors(consB,"B",pcsB,coresets)

    # Angle calculations
    Cvec = as_unit(Bpts.C - Apts.C)
    A1 = as_unit(Apts.V1 - Apts.C)
    A2 = as_unit(Apts.V2 - Apts.C)
    B1 = as_unit(Bpts.V1 - Bpts.C)
    B2 = as_unit(Bpts.V2 - Bpts.C)

    nx = np.cross(A1, Cvec)
    ny = np.cross(Cvec, nx)
    Lp = as_unit([0, np.dot(A1, nx), np.dot(A1, ny)])
    Hp = as_unit([0, np.dot(B1, nx), np.dot(B1, ny)])
    BA = angle(Lp, Hp)
    if np.cross(Lp, Hp)[0] < 0:
        BA = -BA

    BC1 = angle(B1, -Cvec)
    AC1 = angle(A1, Cvec)
    BC2 = angle(B2, -Cvec)
    AC2 = angle(A2, Cvec)
    dc  = np.linalg.norm(Bpts.C - Apts.C)

    out_prefix = Path(out_dir)/Path(input_pdb).stem
    io = PDBIO()
    io.set_structure(input_struct); input_aligned = f"{out_prefix}_aligned_input.pdb"; io.save(input_aligned)
    io.set_structure(consA); consA_out = f"{out_prefix}_consA.pdb"; io.save(consA_out)
    io.set_structure(consB); consB_out = f"{out_prefix}_consB.pdb"; io.save(consB_out)

    generate_pymol_script(input_aligned, consA_out, consB_out, Apts, Bpts, out_prefix)

    return {
        "pdb_name": Path(input_pdb).stem,
        "BA": BA, "BC1": BC1, "AC1": AC1, "BC2": BC2, "AC2": AC2, "dc": dc,
        "input_aligned": input_aligned
    }

# -------------------------------------------------------
# PyMOL visualisation
# -------------------------------------------------------
def generate_pymol_script(input_pdb, consA_pdb, consB_pdb, Apts, Bpts, out_prefix):
    pdb_name = Path(input_pdb).stem
    scale = 15.0
    a1 = Apts.C + scale * (Apts.V1 - Apts.C)
    a2 = Apts.C + scale * (Apts.V2 - Apts.C)
    b1 = Bpts.C + scale * (Bpts.V1 - Bpts.C)
    b2 = Bpts.C + scale * (Bpts.V2 - Bpts.C)

    script = f"""
from pymol import cmd, cgo
cmd.load("{input_pdb}","input_{pdb_name}")
cmd.load("{consA_pdb}","consA_{pdb_name}")
cmd.load("{consB_pdb}","consB_{pdb_name}")
cmd.bg_color("white")
cmd.hide("everything","all")
cmd.show("cartoon","input_{pdb_name}")
cmd.color("blue","input_{pdb_name} and chain A")
cmd.color("green","input_{pdb_name} and chain B")
cmd.set("cartoon_transparency",0.5,"input_{pdb_name}")
cmd.show("cartoon","consA_{pdb_name} or consB_{pdb_name}")
cmd.color("magenta","consA_{pdb_name}")
cmd.color("cyan","consB_{pdb_name}")

cmd.pseudoatom("CAcent_{pdb_name}", pos={list(Apts.C)}, color="red")
cmd.pseudoatom("CBcent_{pdb_name}", pos={list(Bpts.C)}, color="red")
cmd.show("spheres","CAcent_{pdb_name} or CBcent_{pdb_name}")
cmd.set("sphere_scale",1.0,"CAcent_{pdb_name} or CBcent_{pdb_name}")

cmd.pseudoatom("A1end_{pdb_name}", pos={list(a1)})
cmd.pseudoatom("A2end_{pdb_name}", pos={list(a2)})
cmd.pseudoatom("B1end_{pdb_name}", pos={list(b1)})
cmd.pseudoatom("B2end_{pdb_name}", pos={list(b2)})

cmd.load_cgo({add_cgo_arrow(Apts.C,a1,(0.2,0.5,1.0))},"PC1A_{pdb_name}")
cmd.load_cgo({add_cgo_arrow(Apts.C,a2,(0.1,0.8,0.1))},"PC2A_{pdb_name}")
cmd.load_cgo({add_cgo_arrow(Bpts.C,b1,(0.2,0.5,1.0))},"PC1B_{pdb_name}")
cmd.load_cgo({add_cgo_arrow(Bpts.C,b2,(0.1,0.8,0.1))},"PC2B_{pdb_name}")
cmd.load_cgo({add_cgo_arrow(Apts.C,Bpts.C,(0.5,0.0,0.5))},"dc_{pdb_name}")

cmd.angle("AC1_{pdb_name}","A1end_{pdb_name}","CAcent_{pdb_name}","CBcent_{pdb_name}")
cmd.angle("AC2_{pdb_name}","A2end_{pdb_name}","CAcent_{pdb_name}","CBcent_{pdb_name}")
cmd.angle("BC1_{pdb_name}","B1end_{pdb_name}","CBcent_{pdb_name}","CAcent_{pdb_name}")
cmd.angle("BC2_{pdb_name}","B2end_{pdb_name}","CBcent_{pdb_name}","CAcent_{pdb_name}")
cmd.distance("dc_{pdb_name}","CAcent_{pdb_name}","CBcent_{pdb_name}")

cmd.zoom("all")
cmd.png("{out_prefix}_vis.png", dpi=300)
cmd.save("{out_prefix}.pse")
cmd.quit()
"""
    vis_script = f"{out_prefix}_vis.py"
    with open(vis_script, "w") as f:
        f.write(script)
    print(f"✅ PyMOL script written: {vis_script}")

# -------------------------------------------------------
# Main
# -------------------------------------------------------


def values_single_pdb(pdb_file, fv_out, data_path, pcsA, pcsB, coresets, out_dir):
    pdb_name = Path(pdb_file).stem
    fv_path = fv_out/f"{pdb_name}_fv.pdb"
    write_renumbered_fv(str(fv_path), str(pdb_file))
    if not fv_path.exists():
        print(f"Skipping {pdb_name} as Fv file does not exist.")

    result = process(str(fv_path),
                        str(data_path/"consensus_A.pdb"),
                        str(data_path/"consensus_B.pdb"),
                        pcsA, pcsB, coresets, out_dir)
    return result


if __name__ == "__main__":
    data_path = Path("/workspaces/Graphormer/TRangle/data")
    input_folder = Path("/workspaces/Graphormer/TRangle/imgt")
    out_dir = Path("outputs"); out_dir.mkdir(exist_ok=True)
    fv_out = out_dir/"fv"; fv_out.mkdir(exist_ok=True)

    coresets = json.load(open(data_path/"coresets.json"))
    pcsA = np.loadtxt(data_path/"principal_components_alpha.csv")
    pcsB = np.loadtxt(data_path/"principal_components_beta.csv")

    results = []
    result_single=values_single_pdb("/mnt/larry/lilian/DATA/Cory_data/A6/A6prmtop_first_frame.pdb",fv_out, data_path, pcsA, pcsB, coresets, out_dir)
    pdb_name = Path("/workspaces/Graphormer/TRangle/imgt/3vwk.pdb").stem
    print(f"Single PDB result: {result_single}")
    # run PyMOL
    vis_script = Path(out_dir)/f"{pdb_name}_fv_vis.py"
    os.system(f"pymol -cq {vis_script}")

    # clean up temporary pdbs and script
    for suffix in ["_aligned_input.pdb","_consA.pdb","_consB.pdb","_vis.py"]:
        try: os.remove(Path(out_dir)/f"{pdb_name}_fv{suffix}")
        except FileNotFoundError: pass
    quit()

    for pdb_file in input_folder.glob("*.pdb"):
        print(f"Processing {pdb_file.name}...")
        pdb_name = pdb_file.stem
        fv_path = fv_out/f"{pdb_name}_fv.pdb"
        write_renumbered_fv(str(fv_path), str(pdb_file))
        if not fv_path.exists():
            print(f"Skipping {pdb_name} as Fv file does not exist.")
            continue

        result = process(str(fv_path),
                         str(data_path/"consensus_A.pdb"),
                         str(data_path/"consensus_B.pdb"),
                         pcsA, pcsB, coresets, out_dir)
        results.append(result)

        # run PyMOL
        vis_script = Path(out_dir)/f"{pdb_name}_fv_vis.py"
        os.system(f"pymol -cq {vis_script}")

        # clean up temporary pdbs and script
        for suffix in ["_aligned_input.pdb","_consA.pdb","_consB.pdb","_vis.py"]:
            try: os.remove(Path(out_dir)/f"{pdb_name}_fv{suffix}")
            except FileNotFoundError: pass

    # Save CSV
    results_df = pd.DataFrame(results)[["pdb_name","BA","BC1","AC1","BC2","AC2","dc"]]
    results_df.to_csv(out_dir/"angles_results.csv", index=False)
    print(f"✅ Results saved to {out_dir/'angles_results.csv'}")
