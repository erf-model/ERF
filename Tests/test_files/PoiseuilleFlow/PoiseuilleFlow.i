# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =  0    0.  -1.
geometry.prob_hi     =  4.   1.   1.    
amr.n_cell           = 32    4   16

geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"
zhi.type = "NoSlipWall"

# TIME STEP CONTROL
erf.no_substepping     = 1
erf.cfl                = 0.5     # cfl number for hyperbolic system
erf.init_shrink        = 1.0     # scale back initial timestep
erf.change_max         = 1.05    # scale back initial timestep

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = 1000       # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 100        # number of timesteps between plotfiles
erf.plot_vars_1     = density x_velocity y_velocity z_velocity 

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity            = false

erf.les_type         = "None"
erf.molec_diff_type  = "Constant"
erf.rho0_trans       = 1.0
erf.dynamicViscosity = 0.1

erf.use_coriolis = false
erf.abl_driver_type = "PressureGradient"
erf.abl_pressure_grad = -0.2 0. 0.

erf.horiz_spatial_order = 2
erf.vert_spatial_order = 2

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.T_0 = 300.0
prob.prob_type = 10
