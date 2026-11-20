# An anelastic inflow face with no <face>.density: the given theta (301 K, the interior is 300 K)
# is made conserved with the interior base-state density at the same height, which varies by
# 36 % over the column.  After 1000 s the inflow has filled the domain, and theta must be
# 301 K everywhere; a density from another height, or none, puts it far outside.
# This deck: an inflow_outflow face whose velocity and primitive theta come from a
# dirichlet_file; the face must not extrapolate theta from the interior instead.
erf.prob_name = "ABL"
max_step = 100
amrex.fpe_trap_invalid = 1
erf.compressibility_type = Anelastic

# 500 m cells over 4 km: the isentropic base-state density falls from 1.14 to 0.84
geometry.prob_lo     =    0    0    0
geometry.prob_hi     = 4000 1000 4000
amr.n_cell           =    8    2    8
geometry.is_periodic = 0 1 0

xlo.type = "inflow_outflow"
xlo.dirichlet_file = "inflow_file"
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
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta

erf.use_gravity     = true
erf.molec_diff_type = "None"
erf.les_type        = "None"

erf.init_type = Isentropic
prob.rho_0 = 1.2
prob.T_0   = 300.0
prob.U_0   = 10.0
prob.V_0   = 0.0
prob.W_0   = 0.0
