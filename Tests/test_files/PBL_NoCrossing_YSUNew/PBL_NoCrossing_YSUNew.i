# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# PBL height of a column whose bulk Richardson number never reaches Ribcr
# (erf-model/ERF#4458).
#
# A heated surface under a well-mixed 300 K layer that fills the domain, with
# no inversion: the bulk Richardson number stays at or below zero in every
# cell of every column, so the scan reaches the top cell khi without a crossing.  The
# PBL height is then the top of the scan, z(khi) = 31.5 dz = 1008 m, capped at
# 0.9 z(khi) = 907.2 m (WRF bl_ysu.F90 and module_bl_mrf.F, brint = 1).  Before
# the fix, YSUNew gave such a column the floor, max(0.5 dz, 10 m) = 16 m.  The CTest
# entry bounds the pblh of the 2D plotfile.
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

# SURFACE LAYER (heated)
zlo.type = "surface_layer"
erf.most.z0             = 0.1
erf.most.surf_temp_flux = 0.24
erf.most.average_policy = 1
# Keep the PBL height in the surface layer so the 2D plotfile shows it
erf.most.pblh_calc      = "YSU"

zhi.type       = "SlipWall"

# INITIALIZATION: 300 K from the ground to above the domain top, 5 m/s
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_mixed"

# Velocity perturbations only, so the columns differ but theta stays uniform
prob.pert_ref_height = 200.0
prob.pert_deltaU     = 1.0
prob.pert_deltaV     = 1.0
prob.pert_periods_U  = 3.0
prob.pert_periods_V  = 2.0

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
erf.abl_geo_wind    = 5.0 0.0 0.0

# TURBULENCE CLOSURE
erf.les_type = "None"
erf.pbl_type = "YSUNew"
