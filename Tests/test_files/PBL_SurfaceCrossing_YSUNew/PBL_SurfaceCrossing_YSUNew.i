# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# YSUNew PBL height of a column whose surface cell already has a bulk
# Richardson number at or above Ribcr (erf-model/ERF#4458).
#
# A 1 m/s wind over a neutral surface, and a temperature perturbation of
# about 1.5 K in the first cell, of either sign along x (at fixed pressure,
# so theta' is about T' near the ground).  The surface-layer temperature is
# the plane average (erf.most.average_policy = 0), so in the warm columns the
# surface cell's Rib, g (dz/2) theta' / (ws^2 theta), is about 0.8, above
# Ribcr = 0.25: those columns take the floor, max(0.5 dz, 10 m) = 16 m.  In
# the cold columns Rib stays below Ribcr to the top cell: they take
# 0.9 z(khi) = 907.2 m.  Before the fix the interpolation ran with both of
# its ends on the surface cell, 0 / 0, and stopped this run at the
# floating-point trap (amrex.fpe_trap_invalid).  The CTest entry requires the
# pblh of the 2D plotfile to reach both the floor and the top value, so it
# also fails where the trap does not fire and the old code wrote NaN.
erf.prob_name = "ABL"

max_step = 2

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY (one 16x16x32 box: full columns, as the PBL needs)
geometry.prob_extent = 3200  3200  1024
amr.n_cell           =   16    16    32
amr.max_grid_size_x  = 16
amr.max_grid_size_y  = 16
amr.max_grid_size_z  = 32
geometry.is_periodic = 1 1 0

# SURFACE LAYER (neutral, plane-averaged)
zlo.type = "surface_layer"
erf.most.z0             = 0.1
erf.most.surf_temp_flux = 0.0
erf.most.average_policy = 0
# Keep the PBL height in the surface layer so the 2D plotfile shows it
erf.most.pblh_calc      = "YSU"

zhi.type       = "SlipWall"

# INITIALIZATION: 300 K, 1 m/s
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_calm"

# T' = 3 K e^(1/2) (z / 50 m) e^(-(z / 50 m)^2 / 2) cos(2 pi x / Lx): about 1.5 K in the first cell
prob.pert_ref_height = 50.0
prob.pert_deltaT     = 3.0
prob.pert_periods_T  = 1.0

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
erf.plot_vars_1   = density x_velocity y_velocity theta Kmv Khv Lturb
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
erf.abl_geo_wind    = 1.0 0.0 0.0

# TURBULENCE CLOSURE
erf.les_type = "None"
erf.pbl_type = "YSUNew"
