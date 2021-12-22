# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10
stop_time = 900.0

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -12800.   0.    0.
geometry.prob_hi     =  12800. 100. 6400.
amr.n_cell           =  256      4    64     # dx=dy=dz=100 m, Straka et al 1993

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
#erf.fixed_dt       = 1.0     # fixed time step [s] -- from WRF
erf.fixed_dt       = 1.5625e-2  # fixed time step [s] -- Straka et al 1993

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v                = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 1000       # number of timesteps between checkpoints

# PLOTFILES
amr.plot_file       = plt        # root name of plotfile
amr.plot_int        = 3840       # number of timesteps between plotfiles
amr.plot_vars        =  density x_velocity y_velocity z_velocity
amr.derive_plot_vars = pressure theta

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = true
erf.use_coriolis = false
erf.use_rayleigh_damping = false
erf.spatial_order = 2

erf.les_type         = "None"
erf.moledctype       = "Constant"
# diffusion = 75 m^2/s, rho_0 = 1e5/(287*300) = 1.1614401858
erf.dynamicViscosity = 87.108013935 # kg/(m-s)

# PROBLEM PARAMETERS (optional)
prob.T_0 = 300.0
prob.U_0 = 0.0

# SETTING THE TIME STEP
erf.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt
erf.change_max     = 1.05    # multiplier by which dt can change in one time step
erf.init_shrink    = 1.0     # scale back initial timestep

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

## If using the SUNDIALS Submodule, then
## compile with USE_SUNDIALS=TRUE or AMReX_SUNDIALS=ON and
## set strategy here:
#
## integration.sundials.strategy can take the following values:
### NATIVE  = Runge Kutta method controlled by integration.rk
### ERK     = ERKStep from ARKode with SSPRK3 Method
### MRI     = MRIStep from ARKode with Explict Trapezoid Method
### MRITEST = MRIStep from ARKode modified to use no-op inner f0
integration.sundials.strategy = NATIVE
