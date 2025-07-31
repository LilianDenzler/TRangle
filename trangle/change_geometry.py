import math, json
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Superimposer
from pathlib import Path
from SCRIPT_VIS_GEO import write_pseudo_pdb, generate_pymol_script
from collections import namedtuple
from number import write_renumbered_fv
import os
# ========================
# Helpers
# ========================

Points = namedtuple("Points", ["C", "V1", "V2"])

def normalize(v):
    return v / np.linalg.norm(v)

def rotation_matrix(axis, angle_deg):
    axis = normalize(axis)
    theta = math.radians(angle_deg)
    c, s = math.cos(theta), math.sin(theta)
    x, y, z = axis
    C = 1 - c
    return np.array([
        [c+x*x*C, x*y*C-z*s, x*z*C+y*s],
        [y*x*C+z*s, c+y*y*C, y*z*C-x*s],
        [z*x*C-y*s, z*y*C+x*s, c+z*z*C]
    ])

def resolve_plane(a1_deg, a2_deg):
    a1 = math.radians(a1_deg)
    a2 = math.radians(a2_deg)
    v1 = np.array([math.cos(a1), math.sin(a1), 0.0])
    v2 = np.array([math.cos(a2), 0.0, math.sin(a2)])
    return normalize(v1), normalize(v2)


def build_geometry_from_angles(BA, BC1, BC2, AC1, AC2, dc):
    ba_rad = np.radians(BA)
    bc1_rad, bc2_rad = np.radians(BC1), np.radians(BC2)
    ac1_rad, ac2_rad = np.radians(AC1), np.radians(AC2)

    # --- 1. Set up the Foundation ---
    B_centroid = np.array([0.0, 0.0, 0.0])
    A_centroid = np.array([dc, 0.0, 0.0])
    # The normalized inter-centroid vector is the x-axis
    inter_vec_norm = np.array([1.0, 0.0, 0.0])

    # --- 2. Construct B-Domain Axes ---
    # B1 is placed in the xy-plane at angle BC1 to the x-axis
    B1 = np.array([np.cos(bc1_rad), np.sin(bc1_rad), 0.0])

    # B2 must be perpendicular to B1 and at angle BC2 to the x-axis
    # We solve the system of equations: B2.x = cos(BC2) and B2 . B1 = 0
    b2_x = np.cos(bc2_rad)
    # sin(bc1) can't be zero, otherwise B1 is on x-axis and B2 can't be perp to it
    if np.isclose(np.sin(bc1_rad), 0):
        raise ValueError("BC1 cannot be 0 or 180 degrees.")
    b2_y = -b2_x * np.cos(bc1_rad) / np.sin(bc1_rad)

    radicand_b2_z = 1.0 - b2_x**2 - b2_y**2
    if radicand_b2_z < 0:
        raise ValueError("BC1 and BC2 angles are not geometrically compatible.")
    # Choose positive root for the Z component for one solution
    b2_z = np.sqrt(radicand_b2_z)
    B2 = np.array([b2_x, b2_y, b2_z])
    B2 = B2 / np.linalg.norm(B2) # Normalize to be safe

    # --- 3. Construct A-Domain Axes ---
    # A1 is on a cone around the x-axis, rotated by the torsion angle BA
    A1 = np.array([
        -np.cos(ac1_rad),
        np.sin(ac1_rad) * np.cos(ba_rad),
        np.sin(ac1_rad) * np.sin(ba_rad)
    ])

    # A2 construction follows the same logic as B2 construction
    a2_x = -np.cos(ac2_rad)

    # We solve the system A2 . A1 = 0 and |A2|=1 for the other components
    # This is a robust way to solve for the remaining y and z components
    v_perp = np.array([-A1[1], A1[0], 0])
    if np.linalg.norm(v_perp) < 1e-6:
        v_perp = np.array([-A1[2], 0, A1[0]])
    v_perp = v_perp / np.linalg.norm(v_perp)
    v_perp2 = np.cross(A1, v_perp)

    # Solve a 2x2 system for coefficients c1, c2
    # A2 = a2_x*A1 + c1*v_perp + c2*v_perp2 doesn't work.
    # We must solve for y and z directly given A2.x and A2.A1=0
    # A1_x*a2_x + A1_y*a2_y + A1_z*a2_z = 0
    # a2_x^2 + a2_y^2 + a2_z^2 = 1

    # Simplified solution:
    v1 = np.cross(np.array([-1.0, 0, 0]), A1) # Vector in plane perp to A1
    v2 = np.cross(A1, v1)
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)

    # Solve A2 = a*v1_norm + b*v2_norm
    b_coeff = (a2_x - np.dot(np.array([-1.0, 0, 0]), v2_norm) * np.dot(v2_norm, A1)) / np.dot(np.array([-1.0, 0, 0]), v2_norm)
    # Given A1 = (x1, y1, z1) and A2 = (x2, y2, z2) where x2 is known.
    # y1*y2 + z1*z2 = -x1*x2
    # y2^2 + z2^2 = 1 - x2^2
    # This is solving for an intersection of a line and circle in the yz plane.
    y1_p = A1[1]
    z1_p = A1[2]
    c1 = -A1[0] * a2_x # Constant term
    r_sq = 1 - a2_x**2   # Radius squared of the circle

    # Substitute y2 = (c1 - z1_p*z2) / y1_p into circle equation
    if abs(y1_p) > 1e-6:
        # Quadratic equation for z2: (z1_p^2 + y1_p^2) * z2^2 - 2*c1*z1_p*z2 + (c1^2 - r_sq*y1_p^2) = 0
        qa = z1_p**2 + y1_p**2
        qb = -2 * c1 * z1_p
        qc = c1**2 - r_sq * y1_p**2

        radicand_quad = qb**2 - 4 * qa * qc
        if radicand_quad < 0: raise ValueError("AC1/AC2 angles incompatible")

        a2_z = (-qb + np.sqrt(radicand_quad)) / (2 * qa)
        a2_y = (c1 - z1_p * a2_z) / y1_p
    else: # Case where A1 is in xz plane
        a2_z = c1 / z1_p
        radicand_y = r_sq - a2_z**2
        if radicand_y < 0: raise ValueError("AC1/AC2 angles incompatible")
        a2_y = np.sqrt(radicand_y)

    A2 = np.array([a2_x, a2_y, a2_z])
    A2 = A2 / np.linalg.norm(A2)

    return (A_centroid, A1, A2,B_centroid, B1,B2)




def load_structure(path):
    return PDBParser(QUIET=True).get_structure("s", path)

def transform_chain(chain, R, t):
    for atom in chain.get_atoms():
        atom.coord = (R @ atom.coord) + t

def principal_axes(chain):
    cas = np.array([a.coord for a in chain.get_atoms() if a.get_name() == "CA"])
    cas -= cas.mean(axis=0)
    cov = np.cov(cas.T)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    return normalize(vecs[:,order[0]]), normalize(vecs[:,order[1]])

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

def get_coreset_atoms(structure, chain_id, coresets):
    return [
        atom for atom in structure.get_atoms()
        if atom.get_name() == "CA"
        and atom.parent.parent.id == chain_id
        and atom.parent.id[1] in coresets[chain_id]
    ]


# -------------------------------------------------------
# PyMOL visualisation
# -------------------------------------------------------
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


def generate_pymol_script2(input_pdb, cons_pdb, A_C, A_V1, A_V2, B_C, B_V1, B_V2, out_prefix, vis_folder):
    pdb_name = Path(input_pdb).stem
    scale = 10.0
    a1 = A_C + scale * (A_V1 - A_C)
    a2 = A_C + scale * (A_V2 - A_C)
    b1 = B_C + scale * (B_V1 - B_C)
    b2 = B_C + scale * (B_V2 - B_C)

    script = f"""
from pymol import cmd, cgo
cmd.load("{input_pdb}","input_{pdb_name}")
cmd.load("{cons_pdb}","cons_{pdb_name}")
cmd.bg_color("white")
cmd.hide("everything","all")
cmd.show("cartoon","input_{pdb_name}")
cmd.color("blue","input_{pdb_name} and chain A")
cmd.color("green","input_{pdb_name} and chain B")
cmd.set("cartoon_transparency",0.5,"input_{pdb_name}")
cmd.show("cartoon","cons_{pdb_name}")
cmd.color("magenta","cons_{pdb_name}")

cmd.pseudoatom("CAcent_{pdb_name}", pos={list(A_C)}, color="red")
cmd.pseudoatom("CBcent_{pdb_name}", pos={list(B_C)}, color="red")
cmd.show("spheres","CAcent_{pdb_name} or CBcent_{pdb_name}")
cmd.set("sphere_scale",1.0,"CAcent_{pdb_name} or CBcent_{pdb_name}")

cmd.pseudoatom("A1end_{pdb_name}", pos={list(a1)})
cmd.pseudoatom("A2end_{pdb_name}", pos={list(a2)})
cmd.pseudoatom("B1end_{pdb_name}", pos={list(b1)})
cmd.pseudoatom("B2end_{pdb_name}", pos={list(b2)})

cmd.load_cgo({add_cgo_arrow(A_C,a1,(0.2,0.5,1.0))},"PC1A_{pdb_name}")
cmd.load_cgo({add_cgo_arrow(A_C,a2,(0.1,0.8,0.1))},"PC2A_{pdb_name}")
cmd.load_cgo({add_cgo_arrow(B_C,b1,(0.2,0.5,1.0))},"PC1B_{pdb_name}")
cmd.load_cgo({add_cgo_arrow(B_C,b2,(0.1,0.8,0.1))},"PC2B_{pdb_name}")
cmd.load_cgo({add_cgo_arrow(A_C,B_C,(0.5,0.0,0.5))},"dc_{pdb_name}")

cmd.angle("AC1_{pdb_name}","A1end_{pdb_name}","CAcent_{pdb_name}","CBcent_{pdb_name}")
cmd.angle("AC2_{pdb_name}","A2end_{pdb_name}","CAcent_{pdb_name}","CBcent_{pdb_name}")
cmd.angle("BC1_{pdb_name}","B1end_{pdb_name}","CBcent_{pdb_name}","CAcent_{pdb_name}")
cmd.angle("BC2_{pdb_name}","B2end_{pdb_name}","CBcent_{pdb_name}","CAcent_{pdb_name}")
cmd.distance("dc_{pdb_name}","CAcent_{pdb_name}","CBcent_{pdb_name}")

cmd.zoom("all")
cmd.png("{os.path.join(vis_folder,out_prefix+"_vis.png")}", dpi=300)
cmd.save("{os.path.join(vis_folder,out_prefix+"_vis.pse")}")
cmd.quit()
"""
    vis_script = os.path.join(vis_folder,out_prefix+"_vis.py")
    with open(vis_script, "w") as f:
        f.write(script)
    return os.path.join(vis_folder,out_prefix+'_vis.py')


def add_centroid_pseudoatoms(structure, A_C, B_C):
    """
    Inject pseudoatoms (as HETATM) for visualization.
    """
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure
    import numpy as np

    # Create a new structure so we can append easily
    new_struct = Structure("consensus_with_centroids")
    model = Model(0)
    new_struct.add(model)

    # copy over existing chains
    for c in structure[0]:
        model.add(c.copy())

    # helper to create a pseudoatom
    def make_pseudo_atom(name, coord):
        return Atom(
            name, np.array(coord, dtype=float), 1.0, 1.0, ' ', name, 999, 'C'
        )

    # put them in their own chain to avoid confusion
    chainC = Chain("Z")  # new chain for pseudoatoms
    resA = Residue((" ", 900, " "), "CEN", "")
    resA.add(make_pseudo_atom("CA", A_C))
    chainC.add(resA)
    resB = Residue((" ", 901, " "), "CEN", "")
    resB.add(make_pseudo_atom("CB", B_C))
    chainC.add(resB)

    model.add(chainC)
    return new_struct


# -------------------------------------------------------
# main Code
# -------------------------------------------------------
def change_geometry(angles, out_pdb, input_pdb, consA_path, consB_path, coresets, pcsA, pcsB):
    [BA, BC1, BC2, AC1, AC2, dc] = angles

    # Synthetic geometry (target)
    A_C,A_V1, A_V2,B_C, B_V1, B_V2=build_geometry_from_angles(BA, BC1, BC2, AC1, AC2, dc)
    write_pseudo_pdb("geom_points.pdb", A_C,A_V1, A_V2,B_C, B_V1, B_V2)
    generate_pymol_script("geom_points.pdb", A_C,A_V1, A_V2,B_C, B_V1, B_V2)
    os.remove("geom_points.pdb")

    # Load consensus structures
    consA = load_structure(consA_path)
    consB = load_structure(consB_path)
    chainA = [c for c in consA.get_chains() if c.id == "A"][0]
    chainB = [c for c in consB.get_chains() if c.id == "B"][0]

    # ===== Chain A =====
    Apts = centroid_and_vectors(consA, "A", pcsA, coresets)
    A_C_src = Apts.C
    A_V1_src = normalize(Apts.V1 - A_C_src)
    A_V2_src = normalize(Apts.V2 - A_C_src)

    # translate chain A so centroid matches A_C
    tA = A_C - A_C_src
    for atom in chainA.get_atoms():
        atom.coord = atom.coord + tA
    A_V1_t = Apts.V1 + tA
    A_V2_t = Apts.V2 + tA
    A_C_t  = Apts.C  + tA

    # rotate so PCA basis matches synthetic
    R_srcA = np.stack([normalize(A_V1_t - A_C_t),
                       normalize(A_V2_t - A_C_t),
                       np.cross(normalize(A_V1_t - A_C_t),
                                normalize(A_V2_t - A_C_t))], axis=1)
    R_tgtA = np.stack([A_V1, A_V2, np.cross(A_V1, A_V2)], axis=1)
    R_alignA = R_tgtA @ R_srcA.T

    for atom in chainA.get_atoms():
        atom.coord = R_alignA @ (atom.coord - A_C) + A_C
    A_V1_final = R_alignA @ (A_V1_t - A_C) + A_C
    A_V2_final = R_alignA @ (A_V2_t - A_C) + A_C

    # ===== Chain B =====
    Bpts = centroid_and_vectors(consB, "B", pcsB, coresets)
    B_C_src = Bpts.C
    B_V1_src = normalize(Bpts.V1 - B_C_src)
    B_V2_src = normalize(Bpts.V2 - B_C_src)

    # translate chain B so centroid matches B_C
    tB = B_C - B_C_src
    for atom in chainB.get_atoms():
        atom.coord = atom.coord + tB
    B_V1_t = Bpts.V1 + tB
    B_V2_t = Bpts.V2 + tB
    B_C_t  = Bpts.C  + tB

    # rotate so PCA basis matches synthetic
    R_srcB = np.stack([normalize(B_V1_t - B_C_t),
                       normalize(B_V2_t - B_C_t),
                       np.cross(normalize(B_V1_t - B_C_t),
                                normalize(B_V2_t - B_C_t))], axis=1)
    R_tgtB = np.stack([B_V1, B_V2, np.cross(B_V1, B_V2)], axis=1)
    R_alignB = R_tgtB @ R_srcB.T

    for atom in chainB.get_atoms():
        atom.coord = R_alignB @ (atom.coord - B_C) + B_C
    B_V1_final = R_alignB @ (B_V1_t - B_C) + B_C
    B_V2_final = R_alignB @ (B_V2_t - B_C) + B_C

    # ===== Save merged PDB =====
    consA[0].add(chainB)
    # add pseudoatoms for centroids
    new_struct = add_centroid_pseudoatoms(consA, A_C, B_C)
    io = PDBIO()
    io.set_structure(new_struct)
    io.save(out_pdb)
    print("Consensus structure saved to:", out_pdb)
    return new_struct, A_C, A_V1_final, A_V2_final, B_C, B_V1_final, B_V2_final


def apply_transform(structure, rot, tran, chain_name):
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if chain.get_id() != chain_name:
                        continue
                    else:
                        atom.transform(rot, tran)
    return structure


def move_chains_to_geometry(new_consensus, input_pdb, output_pdb, coresets):
    parser = PDBParser(QUIET=True)
    input_struct = parser.get_structure("input", input_pdb)
    # Step 1: superpose input A onto consensus A
    si = Superimposer()
    si.set_atoms(get_coreset_atoms(new_consensus,"A",coresets),
                 get_coreset_atoms(input_struct,"A",coresets))
    rotA, tranA = si.rotran
    input_struct=apply_transform(input_struct, rotA, tranA, "A")

    # Step 2: superpose consensus B onto moved input B
    si2 = Superimposer()
    si2.set_atoms(get_coreset_atoms(new_consensus,"B",coresets),
                  get_coreset_atoms(input_struct,"B",coresets))
    rotB, tranB = si2.rotran
    input_struct=apply_transform(input_struct, rotB, tranB, "B")
    print("Chains moved to new geometry.")
    #save to pdb
    io = PDBIO()
    io.set_structure(input_struct)
    io.save(output_pdb)
    return output_pdb

# ========================
# Main
# ========================
# --- paths ---
data_path = "/workspaces/Graphormer/TRangle/data"
consA_path = f"{data_path}/consensus_A.pdb"
consB_path = f"{data_path}/consensus_B.pdb"
coresets = json.load(open(f"{data_path}/coresets.json"))
pcsA = np.loadtxt(f"{data_path}/principal_components_alpha.csv")
pcsB = np.loadtxt(f"{data_path}/principal_components_beta.csv")

def config(config_file):
    #read config file
    import configparser
    config = configparser.ConfigParser()
    config.read(config_file)
    global BA, BC1, BC2, AC1, AC2, dc
    BA = config.getfloat('Angles', 'BA')
    BC1 = config.getfloat('Angles', 'BC1')
    BC2 = config.getfloat('Angles', 'BC2')
    AC1 = config.getfloat('Angles', 'AC1')
    AC2 = config.getfloat('Angles', 'AC2')
    dc = config.getfloat('Angles', 'dc')
    return BA, BC1, BC2, AC1, AC2, dc

def run(input_pdb, config_file):
    BA, BC1, BC2, AC1, AC2, dc=config(config_file=config_file)
    pdb_name = Path(input_pdb).stem
    out_dir= Path("outputs")
    out_dir.mkdir(exist_ok=True)
    tmp_out= out_dir/pdb_name
    tmp_out.mkdir(exist_ok=True)
    vis_folder = tmp_out / "vis"
    vis_folder.mkdir(exist_ok=True)
    write_renumbered_fv(os.path.join(tmp_out, f"{pdb_name}fv.pdb"), str(input_pdb))
    input_pdb=os.path.join(tmp_out, f"{pdb_name}fv.pdb")
    angles=[BA, BC1, BC2, AC1, AC2, dc]
    new_consensus, A_Cf, A_V1f, A_V2f, B_Cf, B_V1f, B_V2f= change_geometry(angles, os.path.join(tmp_out, f"consensus_oriented.pdb") , input_pdb, consA_path, consB_path, coresets, pcsA, pcsB)
    move_chains_to_geometry(new_consensus, input_pdb, os.path.join(tmp_out, f"{pdb_name}_oriented.pdb"), coresets)
    vis_script=generate_pymol_script2(os.path.join(tmp_out, f"{pdb_name}_oriented.pdb"), os.path.join(tmp_out, f"consensus_oriented.pdb"), A_Cf, A_V1f, A_V2f, B_Cf, B_V1f, B_V2f, pdb_name, vis_folder)

    #run  pymol -cq {vis_script}
    print(f"✅ PyMOL script saved as {vis_script}. Run with:\n   pymol -cq {vis_script}")
    os.system(f"pymol -cq {vis_script}")
    print(f"Output files saved in: {tmp_out}")
    #print angles into txt in the tmp_out folder
    angles_file = tmp_out / "angles.txt"
    with open(angles_file, "w") as f:
        f.write(f"BA: {BA}\nBC1: {BC1}\nBC2: {BC2}\nAC1: {AC1}\nAC2: {AC2}\ndc: {dc}\n")
    print(f"Angles saved to: {angles_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Change TCR geometry based on angles.")
    parser.add_argument('--config', type=str, default='config.ini', help='Path to configuration file.')
    parser.add_argument('--input_pdb', type=str, help='Path to input PDB file.')
    args = parser.parse_args()
    if args.config:
        config_file = args.config
    else:
        config_file = "config.ini"
    if args.input_pdb:
        input_pdb = args.input_pdb
    else:
        raise ValueError("Input PDB file must be specified.")
    run(input_pdb, config_file)


