# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# EB_SquareCylinder_AnelasticAMR with the refined patch narrowed to x, y in [-0.96, 0.96], so its
# lateral coarse/fine boundaries are 0.46 from the cylinder faces.
#
# The level-1 projection is singular (Neumann on every boundary).  On 2 ranks, in the second
# projection of step 1, the BiCGStab bottom solve was asked to reduce a residual of 8e-11 by
# 1e-4, broke down at roundoff, and what it handed back stalled MLMG at resid/bnorm 3.2e-6: the
# run aborted with "MLMG failed".  On 1 rank, with the same boxes, the same breakdown happened
# to be harmless.  solve_with_EB_mlmg now gives the bottom solve an absolute tolerance.  Run on
# 2 ranks (the CTest default) to test this; the CTest entry bounds x_velocity after 4 steps.

erf.prob_name = "EB Square Cylinder"

max_step = 4

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -2.56 -2.56 0.00
geometry.prob_hi     =  5.12  2.56 3.84
amr.n_cell           =  48    32   24
amr.max_grid_size    = 32
amr.blocking_factor  = 4

geometry.is_periodic = 0 0 0

xlo.type = "Inflow"
xhi.type = "Outflow"
xlo.velocity = 1. 0. 0.

ylo.type = "SlipWall"
yhi.type = "SlipWall"
zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.cfl = 0.1

# DIAGNOSTICS & VERBOSITY
erf.v = 1

# REFINEMENT: a patch around the cylinder, below z = 2.56; on 2 ranks its two boxes, split at
# y = 0.32, are on different ranks
amr.max_level = 1
amr.ref_ratio = 2
erf.refinement_indicators = box1
erf.box1.max_level = 1
erf.box1.in_box_lo = -0.96 -0.96 0.00
erf.box1.in_box_hi =  0.96  0.96 2.56

erf.coupling_type = "TwoWay"
erf.cf_width      = 0
erf.cf_set_width  = 0

# EB GEOMETRY: the box reaches below the bottom of the domain
erf.terrain_type   = "EB"
eb2.geometry       = "box"
eb2.small_volfrac  = 1.e-4
eb2.box_lo         = -0.5 -0.5 -1.0
eb2.box_hi         =  0.5  0.5  2.0

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 4
erf.plot_vars_1 = density x_velocity y_velocity z_velocity volfrac theta

# SOLVER CHOICE
erf.anelastic = 1

erf.use_gravity = false

erf.alpha_T = 0.0
erf.alpha_C = 0.0

erf.les_type          = "None"
erf.molec_diff_type   = "Constant"
erf.dynamic_viscosity = 0.004

erf.use_coriolis    = false
erf.abl_driver_type = "None"

erf.init_type = "uniform"
erf.project_initial_velocity = 1

# PROBLEM PARAMETERS
prob.T_0   = 300.0
prob.U_0   = 1.0
prob.rho_0 = 1.0
