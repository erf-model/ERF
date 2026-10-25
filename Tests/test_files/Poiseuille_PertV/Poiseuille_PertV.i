# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Initial state of the Poiseuille channel (prob_type 21, wall-normal direction z) with a
# y-velocity perturbation only: prob.pert_delta_u = 0, prob.pert_delta_v = 2. The
# y-velocity perturbation must appear, and with the two amplitudes swapped it must not
# (Tests/check_v_perturbation.py). Only the initial plotfile is written.
erf.prob_name = "Poiseuille flow"

max_step = 0
amrex.fpe_trap_invalid = 1

erf.anelastic = 1
erf.init_type = uniform

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 6.2831853072  3.1415926536  2
amr.n_cell           = 16  8  16
geometry.is_periodic = 1 1 0

zlo.type = NoSlipWall
zhi.type = NoSlipWall

# TIME STEP CONTROL
erf.fixed_dt = 0.001

# DIAGNOSTICS & VERBOSITY
erf.v = 0
amr.v = 0

amr.max_level = 0

# CHECKPOINT AND PLOT FILES
erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 1
erf.plot_vars_1 = x_velocity y_velocity z_velocity

# SOLVER CHOICE
erf.use_gravity     = false
erf.molec_diff_type = "Constant"
erf.dynamic_viscosity = 0.00253164557

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.T_0   = 300.
prob.prob_type = 21
prob.U_0   = 26.25
prob.pert_delta_u   = 0.0
prob.pert_delta_v   = 2.0
prob.pert_periods_u = 2
prob.pert_periods_v = 2
