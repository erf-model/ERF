# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10
stop_time = 1000.0

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 50. 50. 5000.
amr.n_cell           = 4 4 400

geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt       = 0.5     # fixed time step

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v                = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed
amr.max_grid_size_x = 64      # maximum level number allowed
amr.max_grid_size_y = 64      # maximum level number allowed

# CHECKPOINT FILES
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 100        # number of timesteps between checkpoints

# PLOTFILES
amr.plot_file       = plt        # root name of plotfile
amr.plot_int        = 100      # number of timesteps between plotfiles
amr.plot_vars        =  density x_velocity y_velocity z_velocity
amr.derive_plot_vars = pressure theta

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = false

erf.les_type         = "None"
erf.molec_diff_type  = "Constant"
erf.dynamicViscosity = 5.0

erf.spatial_order = 2

erf.use_coriolis = true
erf.abl_driver_type = "GeostrophicWind"
erf.latitude = 90.
erf.abl_geo_wind = 15.0 0.0 0.0
erf.rotational_time_period = 86164.0900027328

# PROBLEM PARAMETERS (optional)
prob.rho_0 = 1.0
prob.T_0 = 300.0

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
