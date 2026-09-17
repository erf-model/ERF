# Same stationary wall-diffusion state as AnelasticWallDiffusion_X, but run with the
# anelastic midpoint integrator, which leaves the vertical implicit diffusion on.
# The theta walls are in x, so the z boundary conditions stay foextrap and the
# tridiagonal solve is active for theta and for the momenta (fac = 1 0 0: the solve
# happens in the first stage only).  A linear profile is a steady state of the
# diffusion operator, so the implicit solve must reproduce it exactly.
erf.prob_name = "Anelastic Wall Diffusion"
erf.init_type = Uniform
erf.anelastic = 1
erf.use_gravity = false
erf.anelastic_type = MidPoint
erf.molec_diff_type = "ConstantAlpha"
erf.dynamic_viscosity = 0.0
erf.alpha_T = 0.01
erf.alpha_C = 0.0
erf.les_type = "Smagorinsky"
erf.Cs = 0.1
max_step = 2
geometry.prob_lo = 0.0 0.0 0.0
geometry.prob_hi = 2.4 4.8 11.9
amr.n_cell = 8 12 17
geometry.is_periodic = 0 0 0
xlo.type = "NoSlipWall"
xhi.type = "NoSlipWall"
ylo.type = "NoSlipWall"
yhi.type = "NoSlipWall"
zlo.type = "NoSlipWall"
zhi.type = "NoSlipWall"
xlo.theta = 300.0
xhi.theta = 301.0
erf.fixed_dt = 0.001
erf.sum_interval = -1
erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1 = 2
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta
amr.max_level = 0
prob.axis = 0
prob.theta_lo = 300.0
prob.theta_hi = 301.0
