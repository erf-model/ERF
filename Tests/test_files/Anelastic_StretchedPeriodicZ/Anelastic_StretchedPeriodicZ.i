# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# The anelastic projection on a stretched mesh that is periodic in z must stop with a clear
# message (erf-model/ERF#4281).  Only the hybrid FFT solver handles a stretched mesh, and it is
# tridiagonal in z.  The random velocity perturbation gives the initial projection something
# to remove, so the solve cannot be skipped.  The test passes on the abort message.
#
erf.prob_name = "ABL"

max_step = 1
stop_time = 1.0

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = 0. 0. 0.
geometry.prob_hi     = 8. 8. 8.
amr.n_cell           = 8  8  8
geometry.is_periodic = 1 1 1
amr.max_level        = 0

erf.grid_stretching_ratio = 1.05
erf.initial_dz            = 0.8

# TIME INTEGRATION
erf.anelastic = 1
erf.cfl       = 0.5

# DIAGNOSTICS & VERBOSITY
erf.v     = 1
amr.v     = 1

# CHECKPOINT AND PLOT FILES
erf.check_int   = -1
erf.plot_int_1  = -1

# SOLVER CHOICE
erf.use_gravity     = false
erf.les_type        = None
erf.molec_diff_type = None

# INITIAL CONDITIONS
erf.init_type = Uniform
prob.rho_0 = 1.0
prob.T_0   = 300.0
prob.KE_0  = 0.0
prob.A_0   = 0.0
prob.U_0_Pert_Mag = 0.5
prob.V_0_Pert_Mag = 0.5
erf.fix_random_seed = 1
