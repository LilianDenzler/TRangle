
from pymol import cmd, cgo

cmd.load("geom_points.pdb", "geom")
cmd.bg_color("white")
cmd.show("spheres","geom")
cmd.set("sphere_scale",0.5,"geom")

cmd.color("red","geom and name CA1")
cmd.color("blue","geom and name CA2")
cmd.color("blue","geom and name CA3")
cmd.color("green","geom and name CB1")
cmd.color("orange","geom and name CB2")
cmd.color("orange","geom and name CB3")

cmd.load_cgo([
        cgo.CYLINDER,15.200,0.000,0.000,
                     11.517,-5.050,7.806,
                     0.3,
                     0,0,1,
                     0,0,1,
        cgo.CONE,11.517,-5.050,7.806,
                 15.200,0.000,0.000,
                 0.44999999999999996,0.0,
                 0,0,1,
                 0,0,1,1.0
    ],"alpha_V1")
cmd.load_cgo([
        cgo.CYLINDER,15.200,0.000,0.000,
                     24.188,0.213,4.379,
                     0.3,
                     0,0.5,1,
                     0,0.5,1,
        cgo.CONE,24.188,0.213,4.379,
                 15.200,0.000,0.000,
                 0.44999999999999996,0.0,
                 0,0.5,1,
                 0,0.5,1,1.0
    ],"alpha_V2")
cmd.load_cgo([
        cgo.CYLINDER,0.000,0.000,0.000,
                     -3.289,9.444,0.000,
                     0.3,
                     1,0,0,
                     1,0,0,
        cgo.CONE,-3.289,9.444,0.000,
                 0.000,0.000,0.000,
                 0.44999999999999996,0.0,
                 1,0,0,
                 1,0,0,1.0
    ],"beta_V1")
cmd.load_cgo([
        cgo.CYLINDER,0.000,0.000,0.000,
                     4.305,1.499,8.900,
                     0.3,
                     1,0.5,0,
                     1,0.5,0,
        cgo.CONE,4.305,1.499,8.900,
                 0.000,0.000,0.000,
                 0.44999999999999996,0.0,
                 1,0.5,0,
                 1,0.5,0,1.0
    ],"beta_V2")
cmd.load_cgo([
        cgo.CYLINDER,15.200,0.000,0.000,
                     0.000,0.000,0.000,
                     0.3,
                     0,1,0,
                     0,1,0,
        cgo.CONE,0.000,0.000,0.000,
                 15.200,0.000,0.000,
                 0.44999999999999996,0.0,
                 0,1,0,
                 0,1,0,1.0
    ],"C_axis")

cmd.angle("BC1_angle", "geom and name CB2", "geom and name CB1", "geom and name CA1")
cmd.angle("BC2_angle", "geom and name CB3", "geom and name CB1", "geom and name CA1")
cmd.angle("AC1_angle", "geom and name CA2", "geom and name CA1", "geom and name CB1")
cmd.angle("AC2_angle", "geom and name CA3", "geom and name CA1", "geom and name CB1")

cmd.distance("dc_dist","geom and name CA1","geom and name CB1")
cmd.zoom("all")
cmd.png("geometry_visualization.png",dpi=300)
cmd.save("geometry_visualization.pse")
cmd.quit()
