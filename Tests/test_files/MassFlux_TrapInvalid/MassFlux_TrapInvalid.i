# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Anelastic flow in a 2 x 1.5 x 0.8 box (unequal sides, so the x and y areas differ),
# periodic in x and y, with amrex.fpe_trap_invalid = 1.
# MassFlux_Off_TrapInvalid runs one step without a constant mass flux (the default
# erf.const_massflux_u = erf.const_massflux_v = 0): the mass-flux control in
# make_mom_sources then has zero cross-sectional areas, and its guarded divisions must
# not form 0/0 (an optimised build may divide before it tests the area).
# MassFlux_X_Relaxation adds erf.const_massflux_u = 1 over the layer 0.2-0.6 m and checks
# the relaxation of the integrated x mass flux against its closed form (Tests/check_massflux_relaxation.py).
erf.prob_name = "Flow in a Box"

max_step = 1
amrex.fpe_trap_invalid = 1

erf.anelastic = 1
erf.init_type = Uniform

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo = -1. -1.  0.
geometry.prob_hi =  1.  0.5 0.8
amr.n_cell       = 8  8  8
amr.max_grid_size = 8 8 8

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

zlo.theta = 299.
zhi.theta = 280.

# TIME STEP CONTROL
erf.fixed_dt = 0.025

# DIAGNOSTICS & VERBOSITY
erf.v  = 0
amr.v  = 0

# REFINEMENT / REGRIDDING
amr.max_level = 0

# CHECKPOINT AND PLOT FILES
erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = -1

# SOLVER CHOICE
erf.use_gravity     = true
erf.molec_diff_type = "None"
erf.vert_implicit   = false
erf.les_type        = "None"

# PROBLEM PARAMETERS
prob.rho_0        = 1.0
prob.T_0          = 285.
prob.T_0_Pert_Mag = 0.1
