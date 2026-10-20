# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Stable flow over a periodic 2-D ridge on a 3 km mesoscale grid with Smagorinsky2D + MRF
# and no numerical diffusion: the regime where K_h * h^2 can exceed K_v on the slopes.
# Tests/RunSmag2DRidge.cmake runs it with and without the WRF Smagorinsky-2D limits
# (erf.smag2d_slope_limiter, erf.smag2d_kh_cap) and the diffusive time-step check.
erf.prob_name = "ABL"

max_step = 20

amrex.fpe_trap_invalid = 1
# Abort on a NaN on every platform (AMReX traps FPEs only on Linux and macOS)
erf.check_for_nans = 1

# PROBLEM SIZE & GEOMETRY: 32 x 4 cells of 3 km, 40 cells stretched from 50 m (top 12952.8 m)
geometry.prob_extent = 96000  12000  12952.825935499923
amr.n_cell           =    32      4    40
amr.max_grid_size_x  = 16
amr.max_grid_size_y  = 4
amr.max_grid_size_z  = 40
geometry.is_periodic = 1 1 0

erf.grid_stretching_ratio = 1.08
erf.initial_dz            = 50.0

# TERRAIN: Witch-of-Agnesi ridge in x, h(x) = hmax / (1 + (x/L)^2), steepest slope
# 0.65 hmax / L = 0.325 (a 975 m drop across one 3 km cell; alpha = h dx/dz reaches about 19.5
# in the 50 m cells, as the start-up report prints)
erf.terrain_type      = StaticFittedMesh
erf.terrain_smoothing = 2
prob.custom_terrain_type = "WoA"
prob.dir  = 0
prob.hmax = 1500.0
prob.L    = 3000.0

# SURFACE LAYER (adiabatic) AND TOP
zlo.type = "surface_layer"
erf.most.z0 = 0.1
zhi.type = "SlipWall"

# INITIALIZATION: theta + 4 K/km, u from 10 to 20 m/s
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "smag2d_ridge_sounding"

# Absorb the mountain waves in the top 4 km
erf.rayleigh_damp_W   = true
erf.rayleigh_dampcoef = 0.05
erf.rayleigh_zdamp    = 4000.

# TIME STEP CONTROL
erf.fixed_dt           = 10.0
erf.fixed_mri_dt_ratio = 6

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = -1
erf.v            = 1
amr.v            = 0

# REFINEMENT / REGRIDDING
amr.max_level = 0

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 20
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta Kmh Kmv Khv

# SOLVER CHOICE
erf.molec_diff_type = "None"
erf.use_gravity     = true
erf.use_coriolis    = false

erf.dycore_horiz_adv_type = "Upwind_3rd"
erf.dycore_vert_adv_type  = "Upwind_3rd"

# TURBULENCE CLOSURE
erf.les_type = "Smagorinsky2D"
erf.Cs       = 0.25
erf.pbl_type = "MRF"
