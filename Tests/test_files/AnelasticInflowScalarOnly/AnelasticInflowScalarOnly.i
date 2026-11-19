# An anelastic inflow face that gives a scalar but no theta: the scalar enters, made conserved
# with the interior density, and theta is taken from the interior, so it must stay 300 K.  A
# face that gave theta = 0 instead (the missing value of the zero-inflow rule) stopped this run
# with a floating-point error.

erf.prob_name = "ABL"
max_step = 100
amrex.fpe_trap_invalid = 1
erf.compressibility_type = Anelastic

# 500 m cells over 4 km: the isentropic base-state density falls from 1.14 to 0.84
geometry.prob_lo     =    0    0    0
geometry.prob_hi     = 4000 1000 4000
amr.n_cell           =    8    2    8
geometry.is_periodic = 0 1 0

xlo.type = "Inflow"
xlo.velocity = 10.0 0.0 0.0
xlo.scalar = 2.0
xhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.fixed_dt = 10.0
erf.v   = 1
amr.v   = 1
amr.max_level = 0

erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 100
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta scalar

erf.use_gravity     = true
erf.molec_diff_type = "None"
erf.les_type        = "None"

erf.init_type = Isentropic
prob.rho_0 = 1.2
prob.T_0   = 300.0
prob.U_0   = 10.0
prob.V_0   = 0.0
prob.W_0   = 0.0
