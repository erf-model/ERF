# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 20

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  1     1     1
amr.n_cell           = 64     64    4

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.use_lowM_dt    = 1
erf.cfl            = 0.9     # cfl number for hyperbolic system

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
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 20         # number of timesteps between plotfiles
erf.plot_vars_1     = density rhoadv_0 x_velocity y_velocity z_velocity pressure temp theta

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = false

erf.dryscal_horiz_adv_type = "Centered_2nd"
erf.dryscal_vert_adv_type  = "Centered_2nd"
erf.moistscal_horiz_adv_type = "Centered_2nd"
erf.moistscal_vert_adv_type  = "Centered_2nd"

erf.les_type         = "None"
erf.molec_diff_type  = "None"
erf.dynamicViscosity = 0.0

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.T_0   = 1.0
prob.A_0   = 1.0
prob.u_0   = 10.0
prob.v_0   = 5.0
prob.rad_0 = 0.125
prob.uRef  = 0.0
prob.prob_type = 11
