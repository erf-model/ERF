# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# This test demonstrates the use of local timestepping for
# steady-state convergence acceleration in an idealized neutral ABL
# with a single column and vertical grid stretching

max_step = 1000
stop_time = 10000.0  # Run to steady state

amrex.fpe_trap_invalid = 0

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
# Single column: 1x1 in horizontal, 64 cells in vertical
geometry.prob_extent = 10.0 10.0 2000.0
amr.n_cell           = 1 1 64

geometry.is_periodic = 1 1 0

# BOUNDARY CONDITIONS
# Bottom: MOST surface layer
zlo.type = "surface_layer"
erf.most.z0 = 0.1  # roughness length [m]

# Top: Slip wall with fixed theta gradient
zhi.type = "SlipWall"
zhi.theta_grad = 0.003  # 3 K/km lapse rate [K/m]

# VERTICAL GRID STRETCHING
# Use stretched grid to resolve boundary layer
erf.grid_stretching_ratio = 1.03
erf.initial_dz = 10.0
erf.zsurface = 0.0

# TIME STEP CONTROL
erf.fixed_dt = 0.5  # Fixed time step [s]
erf.substepping_type = None

# Enable local timestepping for steady-state convergence
erf.use_local_timestepping = true
erf.smooth_local_dt = true  # Smooth the timestep field for better stability

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 10      # timesteps between computing mass
erf.v            = 1       # verbosity in ERF.cpp
amr.v            = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level = 0  # no refinement

# CHECKPOINT FILES
erf.check_file = chk
erf.check_int  = 500

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 100
erf.plot_vars_1 = density x_velocity y_velocity z_velocity pressure theta

# DATA COLLECTION
erf.data_log = mean_profiles.dat
erf.profile_int = 10

# SOLVER CHOICE
# Use default advection schemes
erf.dycore_horiz_adv_type  = "Centered_2nd"
erf.dycore_vert_adv_type   = "Centered_2nd"
erf.dryscal_horiz_adv_type = "Centered_2nd"
erf.dryscal_vert_adv_type  = "Centered_2nd"

# BACKGROUND CONDITIONS
erf.use_gravity = true
erf.use_coriolis = true
erf.coriolis_3d = false
erf.latitude = 45.0  # degrees

# ABL forcing: geostrophic wind
erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind = 5.0 0.0 0.0  # 5 m/s in x-direction

erf.molec_diff_type = "None"

# TURBULENCE MODELING
# MRF PBL scheme
erf.pbl_type = "MRF"
erf.les_type = "None"

# INITIAL CONDITIONS
erf.init_type = "input_sounding"
erf.sounding_type = "Ideal"

# PROBLEM PARAMETERS
# Small initial perturbations to trigger turbulence
prob.T_0_Pert_Mag = 0.1
prob.pert_rhotheta = false
