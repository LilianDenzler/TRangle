import MDAnalysis as mda
from Bio.PDB import PDBParser, PDBIO
from copy import deepcopy
import os, json, math
from pathlib import Path
from collections import namedtuple
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Superimposer
from .number import write_renumbered_fv, get_renumbered_fv_from_saved
from .new_calc import get_coreset_atoms, centroid_and_vectors, as_unit, angle, apply_transform

def values_single_MD(pdb_file, md_file, out_dir, data_path, pcsA, pcsB, coresets):
    """
    Loop through all frames of an MD trajectory and calculate BA, BC1, AC1, BC2, AC2, dc per frame.
    Returns a pandas DataFrame with columns: frame, BA, BC1, AC1, BC2, AC2, dc
    """
    pdb_name = Path(pdb_file).stem
    fv_path = out_dir / f"{pdb_name}_fv.pdb"
    valid_pair, full_numbering=write_renumbered_fv(str(fv_path), str(pdb_file), fv_only=False)
    if not fv_path.exists():
        print(f"Skipping {pdb_name} as Fv file does not exist.")
        return None

    # Load consensus structures once
    parser = PDBParser(QUIET=True)
    consA = parser.get_structure("consA", data_path / "consensus_A.pdb")
    consB = parser.get_structure("consB", data_path / "consensus_B.pdb")

    # Setup MDAnalysis Universe
    u = mda.Universe(str(pdb_file), str(md_file))

    results = []
    # iterate over trajectory frames
    for i, ts in enumerate(u.trajectory):
        # build a Bio.PDB structure for the current frame
        # Deepcopy the topology structure so we can modify coordinates
        input_struct = parser.get_structure("input", pdb_file)
        for atom_bio, atom_mda in zip(input_struct.get_atoms(), u.atoms):
            atom_bio.coord = atom_mda.position  # overwrite coords with current frame
        input_struct=get_renumbered_fv_from_saved(input_struct, valid_pair, full_numbering)
        #write the modified structure to a temporary file
        #if ts.frame % 100 == 0:
        #    print(f"Processing {pdb_name} frame {i+1}/{len(u.trajectory)}")
        #    pdb_temp = out_dir / f"{pdb_name}_temp_frame_{i}.pdb"
        #    io = PDBIO()
        #    io.set_structure(input_struct)
        #    io.save(str(pdb_temp))
        # Process angles/distances
        result = process_structure(input_struct, consA, consB, pcsA, pcsB, coresets, out_dir, frame=i)
        print(f"Processed frame {i+1}/{len(u.trajectory)}: {result}")
        results.append(result)

    #write results to a csv
    df = pd.DataFrame(results)
    df["frame"] = df.index  # add frame number
    df.to_csv(fv_out / f"{pdb_name}_angles.csv", index=False)
    return df



# A helper wrapper for process() that accepts a structure object to avoid re-parsing each time
def process_structure(input_struct, consA, consB, pcsA, pcsB, coresets, out_dir, frame=None):
    # Deepcopy consensus structures so we can transform them safely
    import copy
    consA_copy = copy.deepcopy(consA)
    consB_copy = copy.deepcopy(consB)
    # Deepcopy input to avoid mutating original
    input_copy = copy.deepcopy(input_struct)

    # Perform alignment and angle calculation
    result = process_inmem(input_copy, consA_copy, consB_copy, pcsA, pcsB, coresets)
    if frame is not None:
        result["frame"] = frame
    return result

# A memory-based variant of your process() (removes file I/O for speed)
def process_inmem(input_struct, consA, consB, pcsA, pcsB, coresets):
    # Step 1: superpose input A onto consensus A
    si = Superimposer()
    si.set_atoms(get_coreset_atoms(consA, "A", coresets),
                 get_coreset_atoms(input_struct, "A", coresets))
    rotA, tranA = si.rotran
    apply_transform(input_struct, rotA, tranA)

    # Step 2: superpose consensus B onto moved input B
    si2 = Superimposer()
    si2.set_atoms(get_coreset_atoms(input_struct, "B", coresets),
                  get_coreset_atoms(consB, "B", coresets))
    rotB, tranB = si2.rotran
    apply_transform(consB, rotB, tranB)

    Apts = centroid_and_vectors(consA, "A", pcsA, coresets)
    Bpts = centroid_and_vectors(consB, "B", pcsB, coresets)

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
    dc = np.linalg.norm(Bpts.C - Apts.C)

    return {"BA": BA, "BC1": BC1, "AC1": AC1, "BC2": BC2, "AC2": AC2, "dc": dc}

# Paths and data loading
data_path = Path("/workspaces/Graphormer/TRangle/data")
coresets = json.load(open(data_path/"coresets.json"))
pcsA = np.loadtxt(data_path/"principal_components_alpha.csv")
pcsB = np.loadtxt(data_path/"principal_components_beta.csv")

def run(pdb_file, md_file, out_dir):
    #make sure out_dir exists
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results=values_single_MD(
        pdb_file=pdb_file,
        md_file=md_file,
        out_dir=out_dir,
        data_path=data_path,
        pcsA=pcsA,
        pcsB=pcsB,
        coresets=coresets
    )
    return results

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calculate angles from MD trajectory")
    parser.add_argument("--pdb_file", type=str, required=True, help="Path to the PDB file")
    parser.add_argument("--md_file", type=str, required=True, help="Path to the MD trajectory file")
    parser.add_argument("--out", type=str, default=str(fv_out), help="Output directory")
    args = parser.parse_args()
    run(pdb_file=args.pdb_file, md_file=args.md_file, out_dir=args.out)
