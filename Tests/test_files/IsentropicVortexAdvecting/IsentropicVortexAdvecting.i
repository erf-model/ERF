# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 1 1 1
geometry.prob_lo     = -12  -12  -1
geometry.prob_hi     =  12   12   1
amr.n_cell           =  48   48   4

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
erf.lo_bc       = "Interior"   "Interior"   "Interior"
erf.hi_bc       = "Interior"   "Interior"   "Interior"

# WHICH PHYSICS
erf.do_hydro = 1

# TIME STEP CONTROL
#erf.cfl            = 0.9     # cfl number for hyperbolic system
#erf.init_shrink    = 1.0     # scale back initial timestep
#erf.change_max     = 1.05    # scale back initial timestep
#erf.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt
erf.fixed_dt        = 0.0005

# DIAGNOSTICS & VERBOSITY
erf.sum_interval    = 1       # timesteps between computing mass
erf.v               = 1       # verbosity in ERF.cpp
amr.v               = 1       # verbosity in Amr.cpp
amr.data_log        = datlog

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
amr.check_int       = 100        # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file       = plt        # root name of plotfile
amr.plot_int        = 10         # number of timesteps between plotfiles
amr.plot_vars        =  density
amr.derive_plot_vars = pressure theta temp x_velocity y_velocity z_velocity

# SOLVER CHOICE
erf.use_state_advection = true
erf.use_momentum_advection = true
erf.use_thermal_diffusion = false
erf.alpha_T = 0.0
erf.use_scalar_diffusion = false
erf.alpha_S = 0.0
erf.use_momentum_diffusion = false
erf.kinematicViscosity = 0.0
erf.use_smagorinsky   = false
erf.use_pressure = true
erf.use_gravity = false
erf.spatial_order = 2

# PROBLEM PARAMETERS
prob.p_inf = 1e5  # reference pressure [Pa]
prob.T_inf = 300. # reference temperature [K]
prob.M_inf = 1.1952286093343936  # freestream Mach number [-]
prob.alpha = 0.7853981633974483  # inflow angle, 0 --> x-aligned [rad]
prob.beta  = 1.1088514254079065 # non-dimensional max perturbation strength [-]
prob.R     = 1.0  # characteristic length scale for grid [m]
prob.sigma = 1.0  # Gaussian standard deviation [-]
#prob.init_periodic = true # initialize a 3x3 array of vortices (8 vortices off-grid)
