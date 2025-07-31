import os
import subprocess
import numpy as np
from Bio.PDB import PDBParser, PDBIO, StructureBuilder
import pathlib
import json
import tempfile

# === CONFIGURATION ===
folder = "/workspaces/Graphormer/TRangle/imgt_variable"
REFERENCE_PDB = "/workspaces/Graphormer/TRangle/imgt_variable/2cdf.pdb"
TALIGN_EXEC = "TMalign"
path = pathlib.Path(__file__).parent
data_path = path.parent / 'data'

with open(data_path / 'coresets.json') as f:
    coresets = json.load(f)

core_alpha = coresets['A']
core_beta = coresets['B']

output_alpha = "consensus_alpha.pdb"
output_beta = "consensus_beta.pdb"

def run_tmalign_chain(mobile_pdb, ref_pdb, chain_id):
    """Align only the given chain and return path to aligned structure with only that chain."""
    with tempfile.NamedTemporaryFile("w+", suffix=".pdb", delete=False) as mob_temp, \
         tempfile.NamedTemporaryFile("w+", suffix=".pdb", delete=False) as ref_temp:

        extract_single_chain(mobile_pdb, chain_id, mob_temp.name)
        extract_single_chain(ref_pdb, chain_id, ref_temp.name)

        subprocess.run(
            [TALIGN_EXEC, mob_temp.name, ref_temp.name, "-o", "TM_sup.pdb"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        os.remove(mob_temp.name)
        os.remove(ref_temp.name)

    if not os.path.exists("TM_sup.pdb"):
        return None

    # Extract only aligned chain A from TM_sup.pdb
    aligned_only = tempfile.NamedTemporaryFile("w+", suffix=".pdb", delete=False)
    with open("TM_sup.pdb") as infile, open(aligned_only.name, "w") as outfile:
        for line in infile:
            if line.startswith("ATOM") and line[21] == 'A':
                outfile.write(line)

    os.remove("TM_sup.pdb")
    os.remove("TM_sup.pdb_all")
    os.remove("TM_sup.pdb_all_atm")
    os.remove("TM_sup.pdb_all_atm_lig")
    os.remove("TM_sup.pdb_atm")
    return aligned_only.name

def extract_single_chain(pdb_path, chain_id, out_path):
    parser = PDBParser(QUIET=True)
    structure = PDBParser(QUIET=True).get_structure("X", pdb_path)
    model = structure[0]

    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found in {pdb_path}")

    chain = model[chain_id]
    io = PDBIO()
    io.set_structure(chain)
    io.save(out_path)

def extract_coords(pdb_path, residue_ids):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_path)
    model = structure[0]
    chain = list(model.get_chains())[0]

    coords = {}
    resid_set = set(residue_ids)

    for res in chain:
        het, resseq, icode = res.id
        if het == ' ' and resseq in resid_set and 'CA' in res:
            coords[resseq] = res['CA'].get_coord()

    return coords

def compute_consensus(folder, chain_id, residue_ids, reference_path):
    pdbs = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(".pdb")]

    # Initialize per-residue list of coordinates
    coord_dict = {resid: [] for resid in residue_ids}

    for pdb in pdbs:
        try:
            aligned_path = run_tmalign_chain(pdb, reference_path, chain_id)
            if aligned_path:
                coords = extract_coords(aligned_path, residue_ids)
                for resid in residue_ids:
                    if resid in coords:
                        coord_dict[resid].append(coords[resid])
                os.remove(aligned_path)
        except Exception as e:
            print(f"Skipping {pdb}: {e}")

    # Compute mean for each residue (only if we have at least one value)
    consensus_coords = {
        resid: np.mean(coord_list, axis=0)
        for resid, coord_list in coord_dict.items()
        if coord_list  # avoid empty lists
    }

    return consensus_coords

def write_ca_only_pdb(coords_dict, chain_id, residue_ids, out_file):
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure("CONS")
    builder.init_model(0)
    builder.init_chain(chain_id)

    for i, res_id in enumerate(residue_ids):
        if res_id not in coords_dict:
            continue  # skip missing
        builder.init_seg(" ")
        builder.init_residue("GLY", " ", res_id, " ")
        builder.init_atom("CA", coords_dict[res_id], 1.0, 1.0, " ", "CA", i, element="C")

    io = PDBIO()
    io.set_structure(builder.get_structure())
    io.save(out_file)

# === MAIN ===
cons_alpha = compute_consensus(folder, "A", core_alpha, REFERENCE_PDB)
cons_beta = compute_consensus(folder, "B", core_beta, REFERENCE_PDB)

write_ca_only_pdb(cons_alpha, "A", core_alpha, output_alpha)
write_ca_only_pdb(cons_beta, "B", core_beta, output_beta)
