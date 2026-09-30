# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# A 100 m Witch-of-Agnesi hill carried by an embedded boundary rather than by
# the mesh, which stays flat.  The grid is 5 m in every direction and the hill
# is centred at x = 300 m, so the terrain is exactly 50 m up at x = 400 m.
#
# Tests/CTestList.cmake uses it for station output over embedded-boundary
# terrain: a station asked for 40 m above the local terrain must sample the
# point asked for at 90 m above z = 0.
#
erf.prob_name = "Flow over Witch of Agnesi hill"

erf.init_type = ConstantDensity

max_step = 100

amrex.fpe_trap_invalid  = 1
amrex.fpe_trap_zero     = 1
amrex.fpe_trap_overflow = 1

fabarray.mfiter_tile_size = 1024 1024 1024

eb2.geometry      = terrain
eb2.small_volfrac = 1.e-4

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =   0.   0.    0.
geometry.prob_hi     = 600.  20.  600.
amr.n_cell           = 120    4   120

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.substepping_type = None
erf.cfl              = 0.5

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = -1
erf.v            = 1
amr.v            = 1

amr.max_level = 0

# CHECKPOINT FILES
erf.check_file = chk
erf.check_int  = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = -1
erf.plot_vars_1 = density x_velocity z_velocity theta

# SOLVER CHOICE
erf.use_gravity  = true
erf.use_coriolis = false
erf.les_type     = "None"

# TERRAIN GRID TYPE
erf.terrain_type      = "EB"
erf.eb_boundary_type  = "NoSlipWall"

erf.molec_diff_type   = "Constant"
erf.dynamic_viscosity = 60.0 # [kg/(m-s)]
erf.alpha_T           = 0.0  # [m^2/s]
erf.theta_ref         = 300.0
erf.eb_diff_constraint_y = true

erf.abl_driver_type   = "PressureGradient"
erf.abl_pressure_grad = -0.02 0. 0.

erf.project_initial_velocity = 0

# PROBLEM PARAMETERS
prob.custom_terrain_type = "WoA"
prob.dir   = 0
prob.T_0   = 300.0
prob.U_0   = 0.0
prob.V_0   = 0.0
prob.rho_0 = 1.16
prob.hmax  = 100.0
prob.L     = 100.0

# STATION OUTPUT on the flank of the hill, where the ground is 50 m up
erf.station_names             = mast
erf.mast.field                = x_velocity theta
erf.mast.x                    = 400.0
erf.mast.y                    = 10.0
erf.mast.height_agl           = 40.0
erf.station_sampling_interval = 1
