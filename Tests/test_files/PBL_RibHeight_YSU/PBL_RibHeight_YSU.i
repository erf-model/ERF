# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# YSU PBL height on a constant-dz mesh, for comparison with the same deck on a
# uniform "stretched" mesh (PBL_RibHeight_YSU_Stretched.i, erf-model/ERF#5167).
#
# A 5 m/s wind over a neutral surface and a stable sounding, theta = 300 K +
# 0.01 K/m z.  The bulk Richardson number is
# Rib(z) = g z (theta_v(z) - theta_v,layer) / (|U|^2 theta_v(z_1)),
# with z the height above the ground; it reaches Ribcr = 0.25 at 146.3 m, a
# few cells up, so the height used in Rib moves the PBL height by several m.
erf.prob_name = "ABL"

max_step = 2

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY (one 8x8x32 box: full columns, as the PBL needs)
geometry.prob_extent = 1600  1600  1024
amr.n_cell           =    8     8    32
amr.max_grid_size_x  = 8
amr.max_grid_size_y  = 8
amr.max_grid_size_z  = 32
geometry.is_periodic = 1 1 0

# SURFACE LAYER (neutral)
zlo.type = "surface_layer"
erf.most.z0             = 0.1
erf.most.surf_temp_flux = 0.0
erf.most.pblh_calc      = "YSU"

zhi.type       = "SlipWall"

# INITIALIZATION: theta = 300 K + 0.01 K/m z, 5 m/s
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_stable"

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
erf.plot_int_1    = 2
erf.plot_vars_1   = density x_velocity y_velocity theta Kmv Khv
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1  = 2
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
