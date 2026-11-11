# A cold, shallow WDM6 moist bubble (BF02 with theta_e = 280 K, q_t = 4 g/kg, 5 km deep): the
# ground is below freezing, so rain, snow and graupel all reach it. After 100 s the largest
# snow_accum is ~0.009 mm. The native C++ WDM6 used to accumulate the fall speeds at the bottom
# level instead of the fallout through it, which gave ~8,600 mm (erf-model/ERF#4364).
# The test requires snow_accum to stay within [0, 1] mm.
erf.prob_name = "Bubble"
erf.init_type = MoistBaseState

max_step  = 200
stop_time = 3600.0

amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_extent = 20000.0 400.0 5000.0
amr.n_cell           = 200     4     50
geometry.is_periodic = 0 1 0
xlo.type = "SlipWall"
xhi.type = "SlipWall"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.fixed_dt = 0.5
erf.fixed_mri_dt_ratio = 4

erf.sum_interval = 1
erf.v = 1
amr.v = 1
amr.max_level = 0

erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 200
erf.plot_vars_1 = density theta temp qv qc qi qrain qsnow qgraup rain_accum snow_accum graupel_accum

erf.use_gravity  = true
erf.use_coriolis = false

erf.dycore_horiz_adv_type    = "Upwind_3rd"
erf.dycore_vert_adv_type     = "Upwind_3rd"
erf.dryscal_horiz_adv_type   = "Upwind_3rd"
erf.dryscal_vert_adv_type    = "Upwind_3rd"
erf.moistscal_horiz_adv_type = "Upwind_3rd"
erf.moistscal_vert_adv_type  = "Upwind_3rd"

erf.les_type       = "None"
erf.pbl_type       = "None"
erf.moisture_model = WDM6
# the native C++ path even in a build with the Fortran bridge (ERF_ENABLE_WDM6_FORT)
erf.use_wdm6_cpp_answer = 1

erf.molec_diff_type   = "ConstantAlpha"
erf.dynamic_viscosity = 0.0
erf.alpha_T           = 0.0
erf.alpha_C           = 0.0

prob.x_c = 10000.0
prob.z_c =  2000.0
prob.x_r =  2000.0
prob.z_r =  2000.0
prob.T_0 =   300.0

prob.do_moist_bubble = true
prob.theta_pert  = 2.0
prob.qt_init     = 0.004
prob.eq_pot_temp = 280.0
