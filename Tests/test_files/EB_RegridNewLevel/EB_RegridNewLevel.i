# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# An EB level made by regrid (erf-model/ERF#5033): flow over a Witch of Agnesi hill of EB
# terrain, at rest at the start and driven by a pressure gradient. amr.max_level = 1, but
# level 1 is tagged where the speed exceeds 1e-8, so it does not exist at the start; the
# regrid after step 1 makes it from level 0 (MakeNewLevelFromCoarse), and a later regrid
# changes its grids (RemakeLevel). The run must take all its steps with the floating-point
# traps on and end with velocities of the size of the flow.
erf.prob_name = "Flow over Witch of Agnesi hill"
max_step = 6
erf.init_type = ConstantDensity

amrex.fpe_trap_invalid=1
amrex.fpe_trap_zero=1
amrex.fpe_trap_overflow=1

fabarray.mfiter_tile_size = 1024 1024 1024

eb2.geometry = terrain
eb2.small_volfrac = 1.e-4

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = 0.   0.   0.
geometry.prob_hi     = 595.0  20.0  600.
amr.n_cell           = 60  4   86
geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.timestepping_type = RK3_NoSubstepping
erf.cfl              = 0.5

# DIAGNOSTICS & VERBOSITY (the test reads the regrid messages from the log)
erf.sum_interval   = 1
erf.v              = 1
amr.v              = 1

# REFINEMENT: level 1 is made by the regrid after step 1
amr.max_level = 1
amr.ref_ratio = 2
erf.regrid_int = 1
erf.refinement_indicators = moving
erf.moving.max_level = 1
erf.moving.field_name = velmag
erf.moving.value_greater = 1.0e-8

erf.coupling_type  = "TwoWay"
erf.cf_width = 0
erf.cf_set_width = 0

# CHECKPOINT AND PLOTFILES
erf.check_int       = -1
erf.plot_file_1     = plt
erf.plot_int_1      = 6
erf.plot_vars_1     = density x_velocity z_velocity theta

# SOLVER CHOICE
erf.use_gravity = true
erf.use_coriolis = false
erf.les_type = "None"

# TERRAIN
erf.terrain_type = "EB"
erf.terrain_smoothing = 0
erf.eb_boundary_type = "NoSlipWall"

# Constant diffusion coefficient
erf.molec_diff_type   = "Constant"
erf.dynamic_viscosity = 60.0 # [kg/(m-s)]
erf.alpha_T           = 0.0 # [m^2/s]
erf.theta_ref = 300.0
erf.eb_diff_constraint_y = true

erf.abl_driver_type = "PressureGradient"
erf.abl_pressure_grad = -0.005 0. 0.

erf.project_initial_velocity = 0

erf.rayleigh_damp_W = true
erf.rayleigh_zdamp = 50.0
erf.rayleigh_dampcoef = 0.25

# PROBLEM PARAMETERS
prob.custom_terrain_type = "WoA"
prob.dir   = 0
prob.T_0   = 300.0
prob.U_0   = 0.0
prob.V_0   = 0.0
prob.rho_0 = 1.16
prob.hmax = 100.0
prob.L = 100.0
