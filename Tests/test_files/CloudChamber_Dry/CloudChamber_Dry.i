# Cloud Chamber Stage 1 dry proof of concept.
# theta and qv below are numerical prescribed state values, not panel
# temperature or a wall-transfer closure.  This deck uses the 2 m x 2 m x 1 m
# reference prism and a deliberately small run for checker-driven CI.
erf.prob_name = "Cloud Chamber"
erf.init_type = Uniform
erf.anelastic = 1
erf.use_gravity = true
erf.vert_implicit = false
erf.molec_diff_type = ConstantAlpha
erf.dynamic_viscosity = 0.0
erf.alpha_T = 0.001
erf.alpha_C = 0.001
erf.les_type = Smagorinsky
erf.Cs = 0.1
erf.fixed_dt = 0.001
erf.sum_interval = -1
erf.check_int = -1
erf.cloud_chamber_budget_interval = 1
erf.plot_file_1 = plt
erf.plot_int_1 = 1
erf.plot_vars_1 = density theta temp pressure x_velocity y_velocity z_velocity

max_step = 2
geometry.prob_lo = 0.0 0.0 0.0
geometry.prob_hi = 2.0 2.0 1.0
geometry.is_periodic = 0 0 0
amr.n_cell = 16 16 8
amr.max_level = 0

prob.theta_bottom = 299.0
prob.theta_top = 280.0
prob.theta_perturbation_amplitude = 0.05
prob.perturbation_mode = deterministic_sine

xlo.type = NoSlipWall
xhi.type = NoSlipWall
ylo.type = NoSlipWall
yhi.type = NoSlipWall
zlo.type = NoSlipWall
zhi.type = NoSlipWall
xlo.theta = 285.0
xhi.theta = 285.0
ylo.theta = 285.0
yhi.theta = 285.0
zlo.theta = 299.0
zhi.theta = 280.0
