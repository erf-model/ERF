# Cloud Chamber Stage 1 SatAdj proof of concept.  qv is a prescribed numerical
# vapor mixing ratio; the deliberately supersaturated value lets the existing
# SatAdj model form qc.  SatAdj does not represent rain, aerosol, or wall
# particle processes.
erf.prob_name = "Cloud Chamber"
erf.init_type = Uniform
erf.anelastic = 1
erf.use_gravity = true
erf.moisture_model = SatAdj
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
erf.plot_vars_1 = density theta temp pressure qv qc x_velocity y_velocity z_velocity

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
prob.qv_bottom = 0.030
prob.qv_top = 0.030

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
xlo.qv = 0.030
xhi.qv = 0.030
ylo.qv = 0.030
yhi.qv = 0.030
zlo.qv = 0.030
zhi.qv = 0.030
