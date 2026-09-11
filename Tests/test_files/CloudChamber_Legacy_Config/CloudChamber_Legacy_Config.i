# Short dry legacy_theta_qv startup compatibility regression.
erf.prob_name = "Cloud Chamber"
erf.init_type = ConstantDensity
erf.anelastic = 1
erf.use_gravity = true
erf.vert_implicit = false
erf.fixed_dt = 0.001
erf.sum_interval = -1
erf.check_int = -1
erf.plot_file_1 = legacy_plt
erf.plot_int_1 = 1
erf.plot_vars_1 = density theta temp pressure x_velocity y_velocity z_velocity

max_step = 0
geometry.prob_lo = 0.0 0.0 0.0
geometry.prob_hi = 2.0 2.0 1.0
geometry.is_periodic = 0 0 0
amr.n_cell = 4 4 4
amr.max_level = 0
amr.max_grid_size = 4

prob.p_inf = 100000.0
prob.T_0 = 292.0
prob.thermodynamic_initialization = legacy_theta_qv
prob.theta_bottom = 300.0
prob.theta_top = 284.0
prob.theta_perturbation_amplitude = 0.02
prob.perturbation_mode = deterministic_sine

xlo.type = NoSlipWall
xlo.theta = 292.0
xhi.type = NoSlipWall
xhi.theta = 292.0
ylo.type = NoSlipWall
ylo.theta = 292.0
yhi.type = NoSlipWall
yhi.theta = 292.0
zlo.type = NoSlipWall
zlo.theta = 300.0
zhi.type = NoSlipWall
zhi.theta = 284.0
