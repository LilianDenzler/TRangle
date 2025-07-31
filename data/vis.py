from pymol import cmd, cgo

cmd.load("aligned.pdb", "structure")
cmd.bg_color('black')
cmd.hide('everything')
cmd.show('cartoon')

L1_vector = [
    cgo.CYLINDER, 125.588, -42.979, -65.407, 133.009, -39.600, -71.196, 0.5,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    cgo.CONE, 133.009, -39.600, -71.196, 125.588, -42.979, -65.407, 0.8, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0
]
H1_vector = [
    cgo.CYLINDER, 127.770, -25.862, -77.133, 120.003, -19.565, -77.262, 0.5,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    cgo.CONE, 120.003, -19.565, -77.262, 127.770, -25.862, -77.133, 0.8, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0
]
centroid_vector = [
    cgo.CYLINDER, 125.588, -42.979, -65.407, 127.770, -25.862, -77.133, 0.5,
    1.0, 0.5, 0.0, 1.0, 0.5, 0.0,
    cgo.CONE, 127.770, -25.862, -77.133, 125.588, -42.979, -65.407, 0.8, 0.0,
    1.0, 0.5, 0.0, 1.0, 0.5, 0.0, 1.0
]

cmd.load_cgo(L1_vector, "L1_vector")
cmd.load_cgo(H1_vector, "H1_vector")
cmd.load_cgo(centroid_vector, "centroid_vector")

cmd.pseudoatom(pos=[125.58790878971877, -42.979076567701796, -65.40731278665172], label="VL_C")
cmd.pseudoatom(pos=[127.7698874180314, -25.86243225537259, -77.13302721351192], label="VH_C")
cmd.pseudoatom(pos=[126.67889810387508, -34.42075441153719, -71.27017000008182], label="HL: 179.63 deg")

cmd.zoom("all")
cmd.ray(1200, 1000)
cmd.png("aligned.png", dpi=300)
cmd.save("visualization.pse")
cmd.quit()
