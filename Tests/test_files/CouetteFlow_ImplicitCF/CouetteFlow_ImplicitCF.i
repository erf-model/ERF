# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Couette flow started from rest, with a refined patch that covers only the lower half of
# the column, so that the diffusion front driven by the lid at zhi crosses a coarse/fine
# boundary in z.
#
# Each column of the fine level ends at that boundary, where the vertical implicit solve
# closes its tridiagonal system with the value on the coarse side.  This test pins the
# answer of that closure: the value must be the coarse data at the stage time, not the
# ghost values left over from the start of the step.  The implicit solve is on in every RK
# stage and the diffusion number is about 6 on the fine level, so the closure carries real
# weight; using the lagged value moves the velocity by about 2e-2.
#
erf.prob_name = "Couette Flow"

max_step = 50
erf.fixed_dt = 0.01          # coarse; the fine level takes two steps of 0.1

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 4  4  16
amr.n_cell           = 4  4  32
amr.max_grid_size_z  = 1024     # columns not split in z, so they end only at the domain or the patch

geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"
zhi.type = "NoSlipWall"

zlo.velocity    = 0.0 0.0 0.0
zhi.velocity    = 2.0 0.0 0.0

# REFINEMENT
amr.max_level = 1
amr.ref_ratio = 2
erf.refinement_indicators = box1
erf.box1.max_level = 1
erf.box1.in_box_lo = 0. 0. 0.
erf.box1.in_box_hi = 4. 4. 8.
erf.coupling_type  = "TwoWay"

# TIME INTEGRATION
erf.substepping_type  = Implicit
erf.vert_implicit_fac = 1        # implicit vertical diffusion in all three RK stages

# The path this test exercises, set explicitly rather than taken from the defaults so that
# a change of default cannot silently move the test off it
erf.implicit_before_substep     = true
erf.implicit_momentum_diffusion = true
erf.implicit_thermal_diffusion  = true

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = -1         # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 50         # number of timesteps between plotfiles
erf.plot_vars_1     = density x_velocity y_velocity z_velocity

# SOLVER CHOICE
erf.use_gravity = false

erf.init_type = Uniform

erf.les_type = None
erf.molec_diff_type   = Constant
erf.dynamic_viscosity = 4.0
erf.alpha_T           = 0.0
erf.alpha_C           = 0.0

# PROBLEM PARAMETERS
prob.prob_type = 1
prob.rho_0 = 1.0
prob.T_0 = 300.0
prob.U_0 = 0.0
prob.V_0 = 0.0
prob.W_0 = 0.0
