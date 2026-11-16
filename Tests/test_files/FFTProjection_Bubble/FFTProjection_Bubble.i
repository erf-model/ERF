# Anelastic warm bubble under an Outflow top, projected with the FFT solvers (issues #5185 and
# #5186).  The tests in CTestList.cmake add grid stretching (hybrid FFT solver) and a refined
# box that reaches only the low end of the periodic x direction.  erf.mg_v = 1 prints the
# divergence before and after every solve, which check_terrain_poisson_divergence.py reads.
geometry.prob_lo     =   0.0   0.0   0.0
geometry.prob_hi     = 320.0 320.0 120.0
amr.n_cell           =  16    16    12
geometry.is_periodic =   1     1     0
amr.max_level        = 0
amrex.fpe_trap_invalid = 1

zlo.type = "SlipWall"
zhi.type = "Outflow"

erf.anelastic = 1
erf.use_fft   = true
erf.mg_v      = 1

erf.molec_diff_type = "None"
erf.les_type        = "None"

erf.prob_name   = "ABL"
erf.init_type   = "Isentropic"   # constant theta = prob.T_0 in hydrostatic balance
erf.use_gravity = true
prob.rho_0 = 1.0
prob.T_0   = 300.0
prob.U_0   = 0.0
prob.V_0   = 0.0
prob.A_0   = 0.0

# A 2 K bubble near the top, so the projection has a correction to make at the top face
prob.bubble_x_c = 160.0
prob.bubble_y_c = 150.0
prob.bubble_z_c = 90.0
prob.bubble_x_r = 80.0
prob.bubble_y_r = 60.0
prob.bubble_z_r = 20.0
prob.bubble_dT  = 2.0

erf.fixed_dt  = 0.05
max_step      = 2
stop_time     = 1000.0

erf.v = 1
amr.v = 0
erf.sum_interval = -1
erf.check_int    = -1
erf.plot_int_1   = -1
