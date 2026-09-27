# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# A 5 m/s westerly over a 100 m hill (Cos4Hill, centred at x = 800 m) on a
# terrain-fitted mesh, with a refined level from the ground up over the hill,
# and station output on the upwind and lee slopes 100 m from the summit, where
# the ground is 53.08 m up.
#
# Tests/CTestList.cmake uses it for
#   - box parity on two levels over terrain (the below-ground ghost cells of a
#     terrain-fitted mesh must agree between boxes, or the refined level, which
#     interpolates from coarse ghost cells, depends on the decomposition), with
#     the refined level from the ground and from 100 m up;
#   - station output over terrain: a station asked for 40 m above the local
#     terrain must sample the point asked for at 93.08 m above z = 0, on this
#     fitted mesh and with the same hill as immersed-forcing terrain.
#
erf.prob_name = "ABL"

max_step  = 20
stop_time = 1.0e6

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =    0.    0.    0.
geometry.prob_hi     = 1600.  800.  600.
amr.n_cell           =   32    16    24

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.terrain_type         = StaticFittedMesh
erf.terrain_smoothing    = 0
prob.custom_terrain_type = "Cos4Hill"
prob.hmax                = 100.0
prob.L                   = 100.0

# TIME STEP CONTROL
# The acoustic substep resolves c dtau / dx = 347 * 0.5/10 / 25 = 0.69 on level 1
erf.fixed_dt           = 0.5
erf.fixed_mri_dt_ratio = 10

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = -1
erf.v              = 1
amr.v              = 1

# REFINEMENT
amr.max_level             = 1
amr.ref_ratio_vect        = 2 2 2
erf.refinement_indicators = box1
erf.box1.max_level        = 1
erf.box1.in_box_lo        =  400.  200.    0.
erf.box1.in_box_hi        = 1200.  600.  300.

# CHECKPOINT FILES
erf.check_file     = chk
erf.check_int      = -1

# PLOTFILES
erf.plot_file_1    = plt
erf.plot_int_1     = 20
erf.plot_vars_1    = density x_velocity y_velocity z_velocity theta

# SOLVER CHOICE
erf.use_gravity     = true
erf.use_coriolis    = false
erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.1
erf.pbl_type        = "None"
erf.abl_driver_type = "None"

# INITIAL CONDITIONS
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "input_sounding"
erf.theta_ref           = 300.0

# STATION OUTPUT on the upwind and lee slopes, 100 m from the summit
erf.station_names             = mast gate
erf.mast.field                = x_velocity y_velocity theta
erf.mast.x                    = 700.0
erf.mast.y                    = 400.0
erf.mast.height_agl           = 40.0
erf.gate.field                = x_velocity z_velocity
erf.gate.x                    = 900.0
erf.gate.y                    = 400.0
erf.gate.height_agl           = 120.0
erf.station_sampling_interval = 1
