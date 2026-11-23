# Implicit Rayleigh damping of w, the default on an anelastic level (erf.rayleigh_damping_type not
# given): anelastic flow over the Schar mountain (U = 10 m/s, N = 0.01 1/s) with Rayleigh damping
# of w at rayleigh_dampcoef * dt = 3.  With erf.rayleigh_damping_type = SlowExplicit the damping of
# w is unstable at the top of the domain there (RK2 amplification 1 - 3 + 9/2 = 2.5 per step) and
# w blows up; damped implicitly it stays at the mountain-wave value, about 1.7 m/s.
erf.prob_name = "Flow over Schar Mountain"
erf.compressibility_type = Anelastic

max_step  = 60
stop_time = 1.0e6
amrex.fpe_trap_invalid = 1

geometry.prob_lo     = -25000.   0.     0.
geometry.prob_hi     =  25000. 500. 21000.
amr.n_cell           =    100    1     70
geometry.is_periodic = 0 1 0
amr.max_level        = 0

xlo.type = "Inflow"
xlo.velocity = 10. 0. 0.
xhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "Outflow"

erf.fixed_dt = 10.0

erf.rayleigh_damp_W   = true
erf.rayleigh_zdamp    = 5000.0
erf.rayleigh_dampcoef = 0.3

erf.init_type = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.sounding_type = Ideal

erf.use_gravity = true
erf.les_type = "None"
erf.molec_diff_type = "None"

erf.terrain_type = StaticFittedMesh
erf.terrain_poisson_solver = mlmg  # needs no FFT
erf.terrain_smoothing = 0
prob.custom_terrain_type = "ScharMountain"
prob.a      = 5000.0
prob.hm     = 250.0
prob.lambda = 4000.0

erf.v = 0
amr.v = 0
erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 60
erf.plot_vars_1 = density x_velocity z_velocity theta
