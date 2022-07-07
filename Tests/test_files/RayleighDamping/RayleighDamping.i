# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 100

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 1 1 0
geometry.prob_extent =  8     8     8
amr.n_cell           = 64     4    64

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
integration.type       = RungeKutta
integration.rk.type    = 3
erf.use_native_mri     = 0
erf.cfl                = 0.9     # cfl number for hyperbolic system
erf.init_shrink        = 1.0     # scale back initial timestep
erf.change_max         = 1.05    # scale back initial timestep

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
erf.plot_file       = plt        # root name of plotfile
erf.plot_int        = 100     # number of timesteps between plotfiles
erf.plot_vars       = x_velocity y_velocity z_velocity
erf.derive_plot_vars = theta

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0

erf.use_gravity = false
erf.use_coriolis = false
erf.use_rayleigh_damping = true

erf.les_type         = "None"
erf.molec_diff_type  = "None"
erf.dynamicViscosity = 0.0

erf.spatial_order = 2

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.T_0   = 1.0
prob.A_0   = 1.0
prob.u_0   = 0.0
prob.v_0   = 0.0
prob.rad_0 = 0.25
prob.z0    = 0.1
prob.zRef  = 80.0
prob.uRef  = 8.0
