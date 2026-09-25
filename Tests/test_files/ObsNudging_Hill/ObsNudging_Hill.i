# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Observation nudging over a hill on a terrain-fitted mesh.
#
# A 5 m/s westerly crosses a 100 m hill centred at x = 800 m.  A met mast on
# the upwind slope and a lidar on the lee slope stand 100 m from the summit,
# where the ground is 53 m up, so a height measured from z = 0 instead of from
# the ground would put their measurements 53 m (five vertical radii) away from
# the cells they should nudge.  The mast measures u, v and theta at 40 m above
# the ground, and the lidar u, v and w with their standard deviations at six
# gates from 40 to 240 m, one of them missing at the first time.  Both files change in time.
# A refined level covers both stations from the ground up, so the nudging and
# the terrain under it are evaluated on two levels.
#
# The regression tests run this deck with nudging and without it, and require
# the station series at the mast and at a lidar gate to end closer to the
# measurements with it; they also check that the answer does not depend on the
# box decomposition or on a restart.
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
erf.mast.height_agl           = 40.0
erf.gate.field                = x_velocity z_velocity
erf.gate.x                    = 900.0
erf.gate.y                    = 400.0
erf.gate.height_agl           = 120.0
erf.station_sampling_interval = 1
