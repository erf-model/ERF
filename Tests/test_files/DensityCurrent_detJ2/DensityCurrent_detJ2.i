# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10
stop_time = 900.0

erf.buoyancy_type = 1

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -12800.   0.    0.
geometry.prob_hi     =  12800. 100. 3200.
amr.n_cell           =  256      4    64     # dx=dy,dz=50m here but z_levels below will make effective dz = 100

geometry.is_periodic = 0 1 0

xlo.type = "Symmetry"
xhi.type = "Outflow"

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt       = 1.0      # fixed time step [s] -- Straka et al 1993
erf.fixed_fast_dt  = 0.25     # fixed time step [s] -- Straka et al 1993

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v                = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = 1000       # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 3840       # number of timesteps between plotfiles
erf.plot_vars_1     = density x_velocity y_velocity z_velocity pressure theta pres_hse dens_hse

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = true
erf.use_coriolis = false
erf.use_rayleigh_damping = false
erf.horiz_spatial_order = 2
erf.vert_spatial_order = 2

erf.les_type         = "None"
erf.molec_diff_type  = "ConstantAlpha"
# diffusion = 75 m^2/s, rho_0 = 1e5/(287*300) = 1.1614401858
erf.dynamicViscosity = 87.108013935 # kg/(m-s)

erf.c_p = 1004.0

# PROBLEM PARAMETERS (optional)
prob.T_0 = 300.0
prob.U_0 = 0.0

# SETTING THE TIME STEP
erf.change_max     = 1.05    # multiplier by which dt can change in one time step
erf.init_shrink    = 1.0     # scale back initial timestep

erf.terrain_z_levels = 0. 100. 200. 300. 400. 500. 600. 700. 800. 900. 1000. 1100. 1200. 1300. 1400. 1500. 1600. 1700. 1800. 1900. 2000. 2100. 2200. 2300. 2400. 2500. 2600. 2700. 2800. 2900. 3000. 3100. 3200. 3300. 3400. 3500. 3600. 3700. 3800. 3900. 4000. 4100. 4200. 4300. 4400. 4500. 4600. 4700. 4800. 4900. 5000. 5100. 5200. 5300. 5400. 5500. 5600. 5700. 5800. 5900. 6000. 6100. 6200. 6300. 6400.

# TERRRAIN GRID TYPE
erf.use_terrain = 1
erf.terrain_smoothing = 1
