import math
import numpy as np
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue, Atom

# -------------------
# Helpers
# -------------------
def normalize(v):
    return v / np.linalg.norm(v)

def rotation_matrix(axis, angle_deg):
    axis = normalize(axis)
    theta = math.radians(angle_deg)
    c, s = math.cos(theta), math.sin(theta)
    x, y, z = axis
    C = 1 - c
    return np.array([
        [c + x*x*C, x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
        [z*x*C - y*s, z*y*C + x*s, c + z*z*C]
    ])

def resolve_plane(a1_deg, a2_deg):
    """Given two bend angles from X-axis, return two perpendicular unit vectors."""
    a1 = math.radians(a1_deg)
    a2 = math.radians(a2_deg)

    # V1: in X-Y plane
    v1 = np.array([math.cos(a1), math.sin(a1), 0.0])

    # V2: also from X axis but rotated out of plane
    # Start in X-Z plane
    v2 = np.array([math.cos(a2), 0.0, math.sin(a2)])
    return normalize(v1), normalize(v2)

def build_geometry(BA, BC1, BC2, AC1, AC2, dc):
    X = np.array([1.0, 0.0, 0.0])
    A_C = np.zeros(3)

    # Alpha
    A_V1, A_V2 = resolve_plane(AC1, AC2)

    # Beta centroid
    B_C = A_C + dc * X

    # Beta: mirror Y and Z to match ABangle convention
    B_V1, B_V2 = resolve_plane(180-BC1, 180-BC2)

    # Apply torsion around X
    Rtwist = rotation_matrix(X, BA)
    B_V1 = normalize(Rtwist @ B_V1)
    B_V2 = normalize(Rtwist @ B_V2)

    return A_C, A_V1, A_V2, B_C, B_V1, B_V2


# -------------------
# PDB writing
# -------------------
def add_atom(chain, name, coord, resseq):
    res = Residue.Residue((' ', resseq, ' '), 'DUM', '')
    atom = Atom.Atom(
        name=name,
        coord=np.array(coord, dtype=float),
        bfactor=1.0,
        occupancy=1.0,
        altloc=' ',
        fullname=f"{name:>4}",
        serial_number=resseq,
        element='C'
    )
    res.add(atom)
    chain.add(res)

def write_pseudo_pdb(filename, A_C, A_V1, A_V2, B_C, B_V1, B_V2, scale=10.0):
    struct = Structure.Structure('geom')
    model = Model.Model(0)
    chain = Chain.Chain('X')
    struct.add(model)
    model.add(chain)

    add_atom(chain, 'CA1', A_C, 1)
    add_atom(chain, 'CA2', A_C + scale*A_V1, 2)
    add_atom(chain, 'CA3', A_C + scale*A_V2, 3)

    add_atom(chain, 'CB1', B_C, 4)
    add_atom(chain, 'CB2', B_C + scale*B_V1, 5)
    add_atom(chain, 'CB3', B_C + scale*B_V2, 6)

    io = PDBIO()
    io.set_structure(struct)
    io.save(filename)

# -------------------
# PyMOL script generation
# -------------------
def cgo_arrow(start, end, color, radius=0.3):
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

def generate_pymol_script(pdb_file, A_C, A_V1, A_V2, B_C, B_V1, B_V2,
                          out_script="vis_geometry.py", scale=10.0):
    a1 = A_C + scale*A_V1
    a2 = A_C + scale*A_V2
    b1 = B_C + scale*B_V1
    b2 = B_C + scale*B_V2

    script = f"""
from pymol import cmd, cgo

cmd.load("{pdb_file}", "geom")
cmd.bg_color("white")
cmd.show("spheres","geom")
cmd.set("sphere_scale",0.5,"geom")

cmd.color("red","geom and name CA1")
cmd.color("blue","geom and name CA2")
cmd.color("blue","geom and name CA3")
cmd.color("green","geom and name CB1")
cmd.color("orange","geom and name CB2")
cmd.color("orange","geom and name CB3")

cmd.load_cgo({cgo_arrow(A_C,a1,(0,0,1))},"alpha_V1")
cmd.load_cgo({cgo_arrow(A_C,a2,(0,0.5,1))},"alpha_V2")
cmd.load_cgo({cgo_arrow(B_C,b1,(1,0,0))},"beta_V1")
cmd.load_cgo({cgo_arrow(B_C,b2,(1,0.5,0))},"beta_V2")
cmd.load_cgo({cgo_arrow(A_C,B_C,(0,1,0))},"C_axis")

cmd.angle("BC1_angle", "geom and name CB2", "geom and name CB1", "geom and name CA1")
cmd.angle("BC2_angle", "geom and name CB3", "geom and name CB1", "geom and name CA1")
cmd.angle("AC1_angle", "geom and name CA2", "geom and name CA1", "geom and name CB1")
cmd.angle("AC2_angle", "geom and name CA3", "geom and name CA1", "geom and name CB1")

cmd.distance("dc_dist","geom and name CA1","geom and name CB1")
cmd.zoom("all")
cmd.png("geometry_visualization.png",dpi=300)
cmd.save("geometry_visualization.pse")
cmd.quit()
"""
    with open(out_script,"w") as f:
        f.write(script)
    print(f"✅ PyMOL script saved as {out_script}. Run with:\n   pymol -cq {out_script}")

# -------------------
# Run example
# -------------------
if __name__=="__main__":
    BA, BC1, BC2, AC1, AC2, dc = 45, 60, 60, 64, 154, 15
    A_C, A_V1, A_V2, B_C, B_V1, B_V2 = build_geometry(BA, BC1, BC2, AC1, AC2, dc)
    write_pseudo_pdb("geom_points.pdb", A_C, A_V1, A_V2, B_C, B_V1, B_V2)
    generate_pymol_script("geom_points.pdb", A_C, A_V1, A_V2, B_C, B_V1, B_V2)
