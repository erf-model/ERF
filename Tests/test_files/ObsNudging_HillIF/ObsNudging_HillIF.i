# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Observation nudging over a hill represented by immersed forcing.
#
# The deck of ObsNudging_Hill on one level, with the same 100 m hill made an
# immersed boundary in a flat mesh instead of the bottom of a terrain-fitted
# one.  The nudging measures the station heights from the terrain surface the
# immersed boundary is built from, so the mast's 40 m and the lidar's 120 m
# gate are 93.08 m and 173.08 m above z = 0 at the stations (the ground there
# is 53.08 m up), and the cells inside the hill are not nudged.  The station
# output below is placed at those absolute heights, since it measures heights
# from the bottom of the mesh.
#
# The regression test runs this deck with nudging and without it and requires
# the series at the mast and at the lidar gate to end closer to the
# measurements with it.
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

erf.terrain_type             = ImmersedForcing
erf.immersed_forcing_substep = true
eb2.small_volfrac            = 0.005
prob.custom_terrain_type = "Cos4Hill"
prob.hmax                = 100.0
prob.L                   = 100.0

# TIME STEP CONTROL
# The acoustic substep resolves c dtau / dx = 347 * 0.5/10 / 50 = 0.35
erf.fixed_dt           = 0.5
erf.fixed_mri_dt_ratio = 10

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = -1
erf.v              = 1
amr.v              = 1

# REFINEMENT
amr.max_level             = 0

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

# OBSERVATION NUDGING
erf.nudging_from_observations     = true
erf.obs_nudging.stations          = mast lidar
erf.obs_nudging.tau               = 5.0
erf.obs_nudging.horizontal_radius = 200.0
erf.obs_nudging.vertical_radius   = 10.0
erf.obs_nudging.sigma_factor      = 1.0
erf.obs_nudging.mast.file         = mast.txt
erf.obs_nudging.mast.x            = 700.0
erf.obs_nudging.mast.y            = 400.0
erf.obs_nudging.lidar.file        = lidar.txt
erf.obs_nudging.lidar.x           = 900.0
erf.obs_nudging.lidar.y           = 400.0

# STATION OUTPUT at the two instruments
erf.station_names             = mast gate
erf.mast.field                = x_velocity y_velocity theta
erf.mast.x                    = 700.0
erf.mast.y                    = 400.0
erf.mast.height_abs           = 93.08
erf.gate.field                = x_velocity z_velocity
erf.gate.x                    = 900.0
erf.gate.y                    = 400.0
erf.gate.height_abs           = 173.08
erf.station_sampling_interval = 1
