# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 20

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  1     1     1    
amr.n_cell           = 16    16    16

geometry.is_periodic = 0 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

xlo.type = "Inflow"
xhi.type = "Outflow"

xlo.velocity = 100. 0. 0.
xlo.density = 1.
xlo.theta = 1.
xlo.scalar = 0.

# TIME STEP CONTROL
erf.use_native_mri = 1
erf.use_lowM_dt = 1
erf.cfl = 0.9

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v                = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = 100        # number of timesteps between checkpoints

# PLOTFILES
erf.plot_int         = 20        # number of timesteps between plotfiles
erf.plot_vars        = density x_velocity y_velocity z_velocity scalar

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 1.0
erf.use_gravity = false

erf.les_type         = "None"
erf.molec_diff_type  = "Constant"
erf.rho0_trans       = 1.0
erf.dynamicViscosity = 0.0

erf.spatial_order = 3

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.A_0 = 1.0
prob.u_0 = 100.0
prob.v_0 = 0.0
prob.uRef  = 0.0

prob.prob_type = 10
