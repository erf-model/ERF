# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 200

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 1 0 1
geometry.prob_lo     =  0     0     0
geometry.prob_hi     = 32     16    4    
amr.n_cell           = 32     16    4

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
erf.lo_bc       = "Interior"   "NoSlipWall"   "Interior"
erf.hi_bc       = "Interior"   "Dirichlet"   "Interior"

yhi.velocity    = 2.0 0.0 0.0  # for Dirichlet BC

# WHICH PHYSICS
erf.do_hydro = 1

# TIME STEP CONTROL
erf.fixed_dt       = 0.1     # fixed time step
#erf.cfl            = 0.9    # cfl number for hyperbolic system
erf.init_shrink    = 1.0     # scale back initial timestep
erf.change_max     = 1.05    # scale back initial timestep
erf.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v                = 1       # verbosity in Amr.cpp
amr.data_log         = datlog

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2 2 2 2 # how often to regrid
amr.blocking_factor = 4       # block factor in grid generation
amr.max_grid_size   = 64
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est

# CHECKPOINT FILES
amr.checkpoint_files_output = 0
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 1000       # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file       = plt        # root name of plotfile
amr.plot_int        = 100        # number of timesteps between plotfiles
amr.plot_vars        = density # u does not get output for some reason
amr.derive_plot_vars = x_velocity y_velocity z_velocity #pressure theta

# SOLVER CHOICE
erf.use_state_advection = true
erf.use_momentum_advection = true
erf.use_thermal_diffusion = false
erf.alpha_T = 0.0
erf.use_scalar_diffusion = false
erf.alpha_C = 1.0
erf.use_momentum_diffusion = true
erf.dynamicViscosity = 0.1
erf.use_smagorinsky   = false
erf.use_pressure = false
erf.use_gravity = false
erf.spatial_order = 2

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.T_0 = 300.0
