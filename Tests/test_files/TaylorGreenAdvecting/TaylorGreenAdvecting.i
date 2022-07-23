# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 1 1 0
geometry.prob_extent = 6.283185307179586476925    6.283185307179586476925    6.283185307179586476925    
amr.n_cell           = 16     16     16

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.use_native_mri     = 1
erf.fixed_dt           = 0.16    # fixed time step
erf.mri_fixed_dt_ratio = 4

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v                = 1       # verbosity in Amr.cpp
amr.data_log         = datlog

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = 100        # number of timesteps between checkpoints

# PLOTFILES
erf.plot_int        = 10        # number of timesteps between plotfiles
erf.plot_vars        =  density x_velocity y_velocity z_velocity pressure temp theta scalar

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = false

erf.les_type         = "None"
erf.molec_diff_type  = "None"
erf.dynamicViscosity = 0.0

erf.spatial_order = 2

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.A_0 = 1.0
prob.V_0 = 1.0
