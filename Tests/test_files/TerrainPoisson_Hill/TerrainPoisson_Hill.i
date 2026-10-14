# The anelastic projection on a terrain-fitted mesh over a 100 m hill, on two levels.
# The refined patch (coarse cells 8-19 in x, 12-23 in y, all of z) has its four
# coarse/fine faces on the flanks of the hill, where the terrain slopes, so the
# refined level is solved with the flux through those faces set to zero
# (erf-model/ERF#4272).  erf.mg_v = 1 prints the divergence after every solve,
# which Tests/check_terrain_poisson_divergence.py reads.
erf.prob_name = "ABL"

max_step  = 6
stop_time = 1.0e6

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -1280. -1280.   0.
geometry.prob_hi     =  1280.  1280. 640.
amr.n_cell           =    32     32   16      # dx = dy = 80 m, dz = 40 m
geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TERRAIN: Witch of Agnesi hill, 100 m high, 300 m half-width, in the middle of the domain
erf.terrain_type         = StaticFittedMesh
erf.terrain_smoothing    = 0
prob.custom_terrain_type = "WoA"
prob.dir                 = 2
prob.hmax                = 100.0
prob.L                   = 300.0

# ANELASTIC; ON A TERRAIN-FITTED MESH THE PROJECTION IS FFT-PRECONDITIONED GMRES.
# erf.use_fft is not set: it is not read on this mesh, and a
# build without FFT stops on it at start-up.
erf.anelastic = 1
erf.mg_v      = 1

# TIME STEP CONTROL
erf.fixed_dt = 2.0

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = -1
erf.v            = 0
amr.v            = 0

# REFINEMENT: z is not refined, and the patch spans all of it
amr.max_level             = 1
amr.ref_ratio_vect        = 2 2 1
erf.refinement_indicators = box1
erf.box1.max_level        = 1
erf.box1.in_box_lo_indices_crse =  8 12  0
erf.box1.in_box_hi_indices_crse = 19 23 15

# CHECKPOINT AND PLOT FILES
erf.check_file  = chk
erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 6
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta

# SOLVER CHOICE
erf.use_gravity     = true
erf.use_coriolis    = false
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.pbl_type        = "None"
erf.abl_driver_type = "None"

# INITIALIZATION: a sheared wind (8 to 12 m/s in x, 3 m/s in y) and a stable theta profile
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "input_sounding"
erf.theta_ref           = 300.0
