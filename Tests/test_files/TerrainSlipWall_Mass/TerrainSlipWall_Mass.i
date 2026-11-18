# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# A closed box (slip walls on all six faces) over terrain that slopes into the x walls: a ridge
# (a Witch of Agnesi in x, 100 m high, half-width 150 m, so about 8 m high and sloping where it
# meets the walls).  A warm bubble rising from above the ridge sets up a flow and a pressure
# perturbation, so nothing enters or leaves and the domain mass ERF prints ("MASS") must not
# change (Tests/RunConservation.cmake).  Fully compressible: the acoustic substeps must keep the
# normal momentum on a wall face at zero, which over sloping terrain they did not do before the
# fix (the terrain term of the horizontal pressure gradient is not zero on the wall face), and
# mass left through the walls.  Tests/CTestList.cmake also runs it with the ridge in y and with
# a radial hill, and without terrain (erf.terrain_type = None, the terrain inputs then unread)
# with theta given on the side walls and with fast Rayleigh damping toward the wind of
# input_sounding_wind (5 / 3 m/s), which the substep without terrain must also keep out of the
# wall faces.
#
erf.prob_name = "ABL"

max_step = 40

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =    0.    0.    0.
geometry.prob_hi     = 1000. 1000.  500.
amr.n_cell           =   32    32    16
amr.max_grid_size    =   16    16    16

geometry.is_periodic = 0 0 0

xlo.type = "SlipWall"
xhi.type = "SlipWall"
ylo.type = "SlipWall"
yhi.type = "SlipWall"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.compressibility_type = FullyCompressible

# TERRAIN
erf.terrain_type         = StaticFittedMesh
prob.custom_terrain_type = WoA
prob.dir                 = 0
prob.hmax                = 100.0
prob.L                   = 150.0

# INITIALIZATION: air at rest, 300 K, and a 10 K bubble
erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"

prob.bubble_dT  = 10.0
prob.bubble_x_c = 400.0
prob.bubble_y_c = 500.0
prob.bubble_z_c = 200.0
prob.bubble_x_r = 150.0
prob.bubble_y_r = 150.0
prob.bubble_z_r = 100.0

# TIME STEP CONTROL: dt = 0.5 s, acoustic substeps of 0.05 s (c dt_fast / dx = 0.55)
erf.fixed_dt           = 0.5
erf.fixed_mri_dt_ratio = 10

# DIAGNOSTICS & VERBOSITY: the integrals at every step, with every digit
erf.sum_interval = 1
erf.v            = 1
amr.v            = 1

amr.max_level = 0

# CHECKPOINT AND PLOT FILES
erf.check_int   = -1
erf.plot_int_1  = -1
erf.plot_vars_1 = density theta

# SOLVER CHOICE
erf.use_gravity   = true
erf.buoyancy_type = 1

erf.les_type        = "None"
erf.molec_diff_type = "None"
