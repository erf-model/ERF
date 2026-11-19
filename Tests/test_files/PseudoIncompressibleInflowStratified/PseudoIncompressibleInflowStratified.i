# A pseudo-incompressible inflow face with no <face>.density, giving theta = 301 K into a 300 K
# domain over a stratified base state (rho0 falls from 1.14 to 0.84 over 4 km).  The ghost
# density is (rho theta)_interior / theta_in, and Rho is typed ext_dir so that the momentum on
# the inflow face uses it: rho theta then enters as rho0 theta0 u and stays rho0 theta0 to the
# tolerance of the solve.  With the interior density on the face it was off by 6e-5 after five
# steps, with the ghost density only in the ghost cells by half that.  The CTests also run it on
# an inflow_outflow face and with the theta from a dirichlet_file (inflow_file).
erf.prob_name = "ABL"
max_step = 5
amrex.fpe_trap_invalid = 1
erf.compressibility_type = PseudoIncompressible

# 500 m cells over 4 km: the isentropic base-state density falls from 1.14 to 0.84
geometry.prob_lo     =    0    0    0
geometry.prob_hi     = 4000 1000 4000
amr.n_cell           =    8    2    8
geometry.is_periodic = 0 1 0

xlo.type = "Inflow"
xlo.velocity = 10.0 0.0 0.0
xlo.theta = 301.0
xhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.fixed_dt = 10.0
erf.v   = 1
erf.sum_interval = 1
amr.v   = 1
amr.max_level = 0

erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = -1
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
