# Active WDM6 cold-wall case. The restart runner checkpoints at step 500 and
# checks the evolved number/mass state before comparing step 1000.
# The test-specific CCN concentration is above WDM6's native 1e8 lower bound
# so evolved Q8 can remain observable across the first resumed step.
# Initial RH is slightly subsaturated so the cold top wall activates only part
# of the column instead of driving every Q8 cell below that native floor.
erf.prob_name = "Cloud Chamber"
erf.init_type = ConstantDensity
erf.anelastic = 1
erf.use_gravity = true
erf.moisture_model = WDM6
erf.vert_implicit = false
erf.molec_diff_type = ConstantAlpha
erf.dynamic_viscosity = 0.0
erf.alpha_T = 0.01
erf.alpha_C = 0.001
erf.les_type = Smagorinsky
erf.Cs = 0.1
erf.fixed_dt = 0.001
erf.sum_interval = -1
erf.check_file = chk
erf.check_int = 500
erf.plot_file_1 = plt
erf.plot_int_1 = 1000
erf.plot_vars_1 = density theta temp pressure qv qc qi qrain qsnow qgraup nc nn nr qt qsat rain_accum snow_accum graup_accum x_velocity y_velocity z_velocity
wdm6.ccn0 = 1000.0e6

max_step = 1000
geometry.prob_lo = 0.0 0.0 0.0
geometry.prob_hi = 2.0 2.0 1.0
geometry.is_periodic = 0 0 0
amr.n_cell = 8 8 8
amr.max_level = 0
amr.max_grid_size = 4

prob.p_inf = 100000.0
prob.T_0 = 284.0
prob.thermodynamic_initialization = physical_temperature_rh
prob.initial_temperature_bottom = 284.0
prob.initial_temperature_top = 284.0
prob.initial_relative_humidity = 0.95
prob.temperature_perturbation_amplitude = 0.0
prob.perturbation_mode = deterministic_sine

xlo.type = NoSlipWall
xhi.type = NoSlipWall
ylo.type = NoSlipWall
yhi.type = NoSlipWall
zlo.type = NoSlipWall
zhi.type = NoSlipWall
xlo.temperature = 284.0
xlo.moisture = dry
xlo.wall_transfer_model = resolved_molecular
xhi.temperature = 284.0
xhi.moisture = dry
xhi.wall_transfer_model = resolved_molecular
ylo.temperature = 284.0
ylo.moisture = dry
ylo.wall_transfer_model = resolved_molecular
yhi.temperature = 284.0
yhi.moisture = dry
yhi.wall_transfer_model = resolved_molecular
zlo.temperature = 284.0
zlo.moisture = dry
zlo.wall_transfer_model = resolved_molecular
zhi.temperature = 270.0
zhi.moisture = dry
zhi.wall_transfer_model = resolved_molecular
