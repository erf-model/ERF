# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 100

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 1 1 1
geometry.prob_extent =  8     8     8
amr.n_cell           = 64     4    64

# TIME STEP CONTROL
erf.cfl            = 0.9     # cfl number for hyperbolic system
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

# CHECKPOINT FILES
amr.checkpoint_files_output = 0
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 100        # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file       = plt        # root name of plotfile
amr.plot_int        = 100     # number of timesteps between plotfiles
amr.plot_vars       = None
amr.derive_plot_vars = theta x_velocity y_velocity z_velocity

# SOLVER CHOICE
erf.use_state_advection = true
erf.use_momentum_advection = true
erf.use_thermal_diffusion = false
erf.alpha_T = 0.0
erf.use_scalar_diffusion = false
erf.alpha_C = 0.0
erf.use_momentum_diffusion = false

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

# INTEGRATION
## integration.type can take on the following values:
## 0 = Forward Euler
## 1 = Explicit Runge Kutta
integration.type = 1

## Explicit Runge-Kutta parameters
#
## integration.rk.type can take the following values:
### 0 = User-specified Butcher Tableau
### 1 = Forward Euler
### 2 = Trapezoid Method
### 3 = SSPRK3 Method
### 4 = RK4 Method
integration.rk.type = 3

## If using a user-specified Butcher Tableau, then
## set nodes, weights, and table entries here:
#
## The Butcher Tableau is read as a flattened,
## lower triangular matrix (but including the diagonal)
## in row major format.
integration.rk.weights = 1
integration.rk.nodes = 0
integration.rk.tableau = 0.0
