# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# YSU cloud PBLH extension with the entrainment term switched off
# (erf-model/ERF#5191).
#
# A heated surface under a 285 K mixed layer, capped by a 10 K inversion whose
# first cell (k = 15, z = 496 m, dz = 32 m) starts supersaturated: the sounding
# gives it 16 g/kg of vapor, and the saturation adjustment of the first step
# leaves about 0.6 g/kg of cloud water there, above erf.ysu_qcloud_threshold
# (0.1 g/kg).  The bulk Richardson number crosses Ribcr in that cell, so it is
# the PBL-top cell, and the PBL height (about 466 m) lies below its top face
# (512 m).
#
# With erf.enable_ysu_cloud_pblh the PBL-top index goes up by one, and the cloudy
# cell takes the in-PBL diffusivity.  Its K profile term vanishes above the PBL
# height (zfac = 1e-8), which leaves
#     Kmv = rho (xkzminm + ckz dz) = rho (0.1 + 0.001 * 32) = 0.132 rho
# against the free-atmosphere floor 0.1 rho without the extension.  The CTest
# entries bound Kmv in that cell.  Before the fix the extension was skipped
# whenever erf.enable_ysu_entrainment = false.
erf.prob_name = "ABL"

max_step = 3

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY (one 8x8x32 box: full columns, as the PBL needs)
geometry.prob_extent = 1600  1600  1024
amr.n_cell           =    8     8    32
amr.max_grid_size_x  = 8
amr.max_grid_size_y  = 8
amr.max_grid_size_z  = 32
geometry.is_periodic = 1 1 0

# SURFACE LAYER (heated)
zlo.type = "surface_layer"
erf.most.z0             = 0.1
erf.most.surf_temp_flux = 0.24
erf.most.average_policy = 1
# Keep the PBL height in the surface layer so the 2D plotfile shows it
erf.most.pblh_calc      = "YSU"

zhi.type       = "SlipWall"

# INITIALIZATION: 285 K mixed layer, inversion with one supersaturated cell, 5 m/s
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_cloud_top"

# TIME STEP CONTROL
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

# PLOTFILES
erf.plot_file_1   = plt
erf.plot_int_1    = 3
erf.plot_vars_1   = density x_velocity y_velocity theta qv qc Kmv Khv Lturb
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1  = 3
erf.plot2d_vars_1 = pblh u_star t_star Olen

# SOLVER CHOICE
erf.molec_diff_type = "None"
erf.use_gravity     = true

erf.use_coriolis           = true
erf.latitude               = 45.0
erf.rotational_time_period = 86455.2516813368

erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind    = 5.0 0.0 0.0

# TURBULENCE CLOSURE
erf.les_type = "None"
erf.pbl_type = "YSU"

# MOISTURE
erf.moisture_model = "SatAdj"

# The two options under test
erf.enable_ysu_entrainment = false
erf.enable_ysu_cloud_pblh  = true
