# Uniform wind turned by Coriolis (inertial oscillation) in a fully periodic box.  The state is
# uniform in space, so every level carries the same values; a regrid at t = 50 s moves the
# fine grids.  A regrid that keeps the time averages leaves them equal on both levels.
max_step = 100
erf.prob_name = "ABL"
erf.init_type = "Uniform"
prob.rho_0 = 1.0
prob.T_0   = 300.0
prob.U_0   = 10.0
prob.V_0   = 0.0

geometry.prob_extent = 1000. 1000. 1000.
geometry.is_periodic = 1 1 1
amr.n_cell           = 16 16 16
amr.max_level        = 1
amr.ref_ratio        = 2
erf.dt_ref_ratio     = 1
erf.regrid_int       = 1

erf.refinement_indicators = box1 box2
erf.box1.in_box_lo = 250. 250. 250.
erf.box1.in_box_hi = 500. 500. 500.
erf.box1.end_time  = 50.
erf.box2.in_box_lo = 375. 375. 375.
erf.box2.in_box_hi = 750. 750. 750.
erf.box2.start_time = 50.

erf.use_gravity  = false
erf.use_coriolis = true
erf.latitude     = 90.
erf.rotational_time_period = 600.
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.pbl_type        = "None"
erf.fixed_dt        = 1.0

erf.time_avg_vel      = true
erf.compute_mean_vars = true

erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 100
erf.plot_vars_1 = x_velocity y_velocity u_t_avg v_t_avg umag_t_avg u_mean v_mean uu_mean uw_mean theta_mean
erf.v = 0
