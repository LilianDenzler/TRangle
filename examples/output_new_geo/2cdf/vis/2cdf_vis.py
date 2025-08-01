
from pymol import cmd, cgo
cmd.load("output_new_geo/2cdf/2cdf_oriented.pdb","input_2cdf_oriented")
cmd.load("output_new_geo/2cdf/consensus_oriented.pdb","cons_2cdf_oriented")
cmd.bg_color("white")
cmd.hide("everything","all")
cmd.show("cartoon","input_2cdf_oriented")
cmd.color("blue","input_2cdf_oriented and chain A")
cmd.color("green","input_2cdf_oriented and chain B")
cmd.set("cartoon_transparency",0.5,"input_2cdf_oriented")
cmd.show("cartoon","cons_2cdf_oriented")
cmd.color("magenta","cons_2cdf_oriented")

cmd.pseudoatom("CAcent_2cdf_oriented", pos=[15.2, 0.0, 0.0], color="red")
cmd.pseudoatom("CBcent_2cdf_oriented", pos=[0.0, 0.0, 0.0], color="red")
cmd.show("spheres","CAcent_2cdf_oriented or CBcent_2cdf_oriented")
cmd.set("sphere_scale",1.0,"CAcent_2cdf_oriented or CBcent_2cdf_oriented")

cmd.pseudoatom("A1end_2cdf_oriented", pos=[11.517131733692704, -5.049959295100113, 7.806048574662422])
cmd.pseudoatom("A2end_2cdf_oriented", pos=[24.187940573750208, 0.21336915651470284, 4.378515771449429])
cmd.pseudoatom("B1end_2cdf_oriented", pos=[-3.2886664603581632, 9.443763763661277, 5.230899662991048e-08])
cmd.pseudoatom("B2end_2cdf_oriented", pos=[4.30511084761935, 1.4991982769639511, 8.900473044600139])

cmd.load_cgo([
        cgo.CYLINDER,15.200,0.000,0.000,
                     11.517,-5.050,7.806,
                     0.3,
                     0.2,0.5,1.0,
                     0.2,0.5,1.0,
        cgo.CONE,11.517,-5.050,7.806,
                 15.200,0.000,0.000,
                 0.44999999999999996,0.0,
                 0.2,0.5,1.0,
                 0.2,0.5,1.0,1.0
    ],"PC1A_2cdf_oriented")
cmd.load_cgo([
        cgo.CYLINDER,15.200,0.000,0.000,
                     24.188,0.213,4.379,
                     0.3,
                     0.1,0.8,0.1,
                     0.1,0.8,0.1,
        cgo.CONE,24.188,0.213,4.379,
                 15.200,0.000,0.000,
                 0.44999999999999996,0.0,
                 0.1,0.8,0.1,
                 0.1,0.8,0.1,1.0
    ],"PC2A_2cdf_oriented")
cmd.load_cgo([
        cgo.CYLINDER,0.000,0.000,0.000,
                     -3.289,9.444,0.000,
                     0.3,
                     0.2,0.5,1.0,
                     0.2,0.5,1.0,
        cgo.CONE,-3.289,9.444,0.000,
                 0.000,0.000,0.000,
                 0.44999999999999996,0.0,
                 0.2,0.5,1.0,
                 0.2,0.5,1.0,1.0
    ],"PC1B_2cdf_oriented")
cmd.load_cgo([
        cgo.CYLINDER,0.000,0.000,0.000,
                     4.305,1.499,8.900,
                     0.3,
                     0.1,0.8,0.1,
                     0.1,0.8,0.1,
        cgo.CONE,4.305,1.499,8.900,
                 0.000,0.000,0.000,
                 0.44999999999999996,0.0,
                 0.1,0.8,0.1,
                 0.1,0.8,0.1,1.0
    ],"PC2B_2cdf_oriented")
cmd.load_cgo([
        cgo.CYLINDER,15.200,0.000,0.000,
                     0.000,0.000,0.000,
                     0.3,
                     0.5,0.0,0.5,
                     0.5,0.0,0.5,
        cgo.CONE,0.000,0.000,0.000,
                 15.200,0.000,0.000,
                 0.44999999999999996,0.0,
                 0.5,0.0,0.5,
                 0.5,0.0,0.5,1.0
    ],"dc_2cdf_oriented")

cmd.angle("AC1_2cdf_oriented","A1end_2cdf_oriented","CAcent_2cdf_oriented","CBcent_2cdf_oriented")
cmd.angle("AC2_2cdf_oriented","A2end_2cdf_oriented","CAcent_2cdf_oriented","CBcent_2cdf_oriented")
cmd.angle("BC1_2cdf_oriented","B1end_2cdf_oriented","CBcent_2cdf_oriented","CAcent_2cdf_oriented")
cmd.angle("BC2_2cdf_oriented","B2end_2cdf_oriented","CBcent_2cdf_oriented","CAcent_2cdf_oriented")
cmd.distance("dc_2cdf_oriented","CAcent_2cdf_oriented","CBcent_2cdf_oriented")

cmd.zoom("all")
cmd.png("output_new_geo/2cdf/vis/2cdf_vis.png", dpi=300)
cmd.save("output_new_geo/2cdf/vis/2cdf_vis.pse")
cmd.quit()
