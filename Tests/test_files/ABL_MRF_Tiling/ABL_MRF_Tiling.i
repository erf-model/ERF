# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# MRF on boxes wider than the MFIter tile size.  Tests/RunTilingParity.cmake
# runs this deck with fabarray.mfiter_tile_size = 1024000 8 8 (the AMReX CPU
# default, so each 16-cell-wide box is split into two tiles in y) and again with
# tiling off, and requires identical plotfiles.  Do not set mfiter_tile_size
# here: the script supplies it.  The same deck runs YSUNew, and legacy YSU with
# a cooled surface (it aborts in unstable conditions), through RUNTIME_OPTIONS.
erf.prob_name = "ABL"

max_step = 10

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 6400  6400  1024
amr.n_cell           =   32    32    32
# Four 16x16x32 boxes: full columns (MRF requires it), wider than a tile in y
amr.max_grid_size_x  = 16
amr.max_grid_size_y  = 16
amr.max_grid_size_z  = 32
geometry.is_periodic = 1 1 0

# SURFACE LAYER (heated, so the countergradient terms are active)
zlo.type = "surface_layer"
erf.most.z0             = 0.1
erf.most.surf_temp_flux = 0.24
# Keep the MRF PBL height in the surface layer so the 2D plotfile shows it
erf.most.pblh_calc      = "MRF"

zhi.type       = "SlipWall"
zhi.theta_grad = 0.003

# INITIALIZATION
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "mrf_sounding_unstable"

# Deterministic (position-only) perturbations so the columns differ and the
# initial state does not depend on the decomposition
prob.pert_ref_height = 200.0
prob.pert_deltaU     = 1.0
prob.pert_deltaV     = 1.0
prob.pert_periods_U  = 3.0
prob.pert_periods_V  = 2.0
prob.pert_deltaT     = 0.5
prob.pert_periods_T  = 2.0

# TIME STEP CONTROL (dx = 200 m: at dx = 100 m with ratio 4 the acoustic
# substep is unstable and the run goes to NaN by step 5)
erf.fixed_dt           = 1.0
erf.fixed_mri_dt_ratio = 6

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.v            = 1
amr.v            = 1

# REFINEMENT / REGRIDDING
amr.max_level = 0

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES (file prefixes are set by the script)
erf.plot_file_1   = plt
erf.plot_int_1    = 10
erf.plot_vars_1   = density x_velocity y_velocity z_velocity pressure theta Kmv Khv Lturb
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1  = 10
erf.plot2d_vars_1 = pblh u_star t_star Olen

# SOLVER CHOICE
erf.molec_diff_type = "None"
erf.use_gravity     = true

erf.use_coriolis           = true
erf.latitude               = 45.0
erf.rotational_time_period = 86455.2516813368

erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind    = 15.0 0.0 0.0

# TURBULENCE CLOSURE
erf.les_type = "None"
erf.pbl_type = "MRF"
