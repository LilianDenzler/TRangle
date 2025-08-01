
from pymol import cmd, cgo
cmd.load("/workspaces/Graphormer/TRangle/examples/output/vis/aligned_input.pdb","input_aligned_input")
cmd.load("/workspaces/Graphormer/TRangle/examples/output/vis/consA.pdb","consA_aligned_input")
cmd.load("/workspaces/Graphormer/TRangle/examples/output/vis/consB.pdb","consB_aligned_input")
cmd.bg_color("white")
cmd.hide("everything","all")
cmd.show("cartoon","input_aligned_input")
cmd.color("blue","input_aligned_input and chain A")
cmd.color("green","input_aligned_input and chain B")
cmd.set("cartoon_transparency",0.5,"input_aligned_input")
cmd.show("cartoon","consA_aligned_input or consB_aligned_input")
cmd.color("magenta","consA_aligned_input")
cmd.color("cyan","consB_aligned_input")

cmd.pseudoatom("CAcent_aligned_input", pos=[110.54427337646486, -9.201766967773432, -72.74709320068361], color="red")
cmd.pseudoatom("CBcent_aligned_input", pos=[101.04493870127658, -20.77018125168697, -74.77665571830558], color="red")
cmd.show("spheres","CAcent_aligned_input or CBcent_aligned_input")
cmd.set("sphere_scale",1.0,"CAcent_aligned_input or CBcent_aligned_input")

cmd.pseudoatom("A1end_aligned_input", pos=[117.4643999040127, -20.710294365882866, -66.06385669112207])
cmd.pseudoatom("A2end_aligned_input", pos=[120.33423614501955, 0.30300724506378707, -66.51685228943826])
cmd.pseudoatom("B1end_aligned_input", pos=[86.54932542431811, -20.957906811495285, -70.92412961981582])
cmd.pseudoatom("B2end_aligned_input", pos=[103.10111208546618, -8.470423651856343, -66.44070325753974])

cmd.load_cgo([
        cgo.CYLINDER,110.544,-9.202,-72.747,
                     117.464,-20.710,-66.064,
                     0.3,
                     0.2,0.5,1.0,
                     0.2,0.5,1.0,
        cgo.CONE,117.464,-20.710,-66.064,
                 110.544,-9.202,-72.747,
                 0.44999999999999996,0.0,
                 0.2,0.5,1.0,
                 0.2,0.5,1.0,1.0
    ],"PC1A_aligned_input")
cmd.load_cgo([
        cgo.CYLINDER,110.544,-9.202,-72.747,
                     120.334,0.303,-66.517,
                     0.3,
                     0.1,0.8,0.1,
                     0.1,0.8,0.1,
        cgo.CONE,120.334,0.303,-66.517,
                 110.544,-9.202,-72.747,
                 0.44999999999999996,0.0,
                 0.1,0.8,0.1,
                 0.1,0.8,0.1,1.0
    ],"PC2A_aligned_input")
cmd.load_cgo([
        cgo.CYLINDER,101.045,-20.770,-74.777,
                     86.549,-20.958,-70.924,
                     0.3,
                     0.2,0.5,1.0,
                     0.2,0.5,1.0,
        cgo.CONE,86.549,-20.958,-70.924,
                 101.045,-20.770,-74.777,
                 0.44999999999999996,0.0,
                 0.2,0.5,1.0,
                 0.2,0.5,1.0,1.0
    ],"PC1B_aligned_input")
cmd.load_cgo([
        cgo.CYLINDER,101.045,-20.770,-74.777,
                     103.101,-8.470,-66.441,
                     0.3,
                     0.1,0.8,0.1,
                     0.1,0.8,0.1,
        cgo.CONE,103.101,-8.470,-66.441,
                 101.045,-20.770,-74.777,
                 0.44999999999999996,0.0,
                 0.1,0.8,0.1,
                 0.1,0.8,0.1,1.0
    ],"PC2B_aligned_input")
cmd.load_cgo([
        cgo.CYLINDER,110.544,-9.202,-72.747,
                     101.045,-20.770,-74.777,
                     0.3,
                     0.5,0.0,0.5,
                     0.5,0.0,0.5,
        cgo.CONE,101.045,-20.770,-74.777,
                 110.544,-9.202,-72.747,
                 0.44999999999999996,0.0,
                 0.5,0.0,0.5,
                 0.5,0.0,0.5,1.0
    ],"dc_aligned_input")

cmd.angle("AC1_aligned_input","A1end_aligned_input","CAcent_aligned_input","CBcent_aligned_input")
cmd.angle("AC2_aligned_input","A2end_aligned_input","CAcent_aligned_input","CBcent_aligned_input")
cmd.angle("BC1_aligned_input","B1end_aligned_input","CBcent_aligned_input","CAcent_aligned_input")
cmd.angle("BC2_aligned_input","B2end_aligned_input","CBcent_aligned_input","CAcent_aligned_input")
cmd.distance("dc_aligned_input","CAcent_aligned_input","CBcent_aligned_input")

cmd.zoom("all")
cmd.png("output/vis/vis.png", dpi=300)
cmd.save("output/vis/vis.pse")
cmd.quit()
