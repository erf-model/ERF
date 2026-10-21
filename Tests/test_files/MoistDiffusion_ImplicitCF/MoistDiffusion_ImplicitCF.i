# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# A moist layer (qv = 10 g/kg above z = 12 m, dry below) diffusing down through a
# coarse/fine boundary in z: the refined patch covers only the lower half of the column.
#
# The vertical implicit solve for moisture (ERF_ImplicitPost.H) runs on the state whose
# valid region has been advanced to the stage time, and closes each fine column at the
# coarse/fine boundary with the value on the coarse side.  This test pins the answer of
# that closure: the value must be the coarse data at the stage time, not the ghost values
# left over from the previous stage.  The companion test CouetteFlow_ImplicitCF covers the
# solve for theta and the momenta in slow_rhs_fun_pre.  Using the lagged value moves qv
# by about 1e-5 here.
#
erf.prob_name = "ABL"

max_step = 200
erf.fixed_dt = 0.005         # coarse; the fine level takes two steps of 0.025

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 4  4  16
amr.n_cell           = 4  4  32
amr.max_grid_size_z  = 1024     # columns not split in z, so they end only at the domain or the patch

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

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

# The path this test exercises, set explicitly rather than taken from the defaults
erf.implicit_before_substep     = true
erf.implicit_moisture_diffusion = true

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1
erf.v              = 1
amr.v              = 1

# CHECKPOINT FILES
erf.check_file      = chk
erf.check_int       = -1

# PLOTFILES
erf.plot_file_1     = plt
erf.plot_int_1      = 200
erf.plot_vars_1     = density qv

# SOLVER CHOICE
erf.use_gravity = false

erf.init_type = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.sounding_type = ConstantDensity

erf.moisture_model = "MoistNoCondensation"

erf.les_type = None
erf.molec_diff_type   = Constant
erf.dynamic_viscosity = 0.0
erf.alpha_T           = 0.0
erf.alpha_C           = 4.0

# PROBLEM PARAMETERS
prob.KE_0 = 0.0
# No passive scalar: its diffusion (alpha_C) is explicit, and at this time step it would
# be unstable on the fine level (diffusion number 1.6)
prob.A_0 = 0.0
