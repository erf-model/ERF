# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# A scalar blob at rest diffusing through the faces of a refined patch, for the conservation tests
# of the two-way reflux.  The patch covers the middle half of the domain in x and in y and the
# lower half in z, so it has coarse/fine faces normal to all three directions, and two-way
# coupling refluxes the scalar across them.  With slip walls and periodic sides nothing leaves the
# domain, so the composite integral of rho*scalar that ERF prints ("RHO SCALAR SL/ML") must not
# change: that holds only if the diffusive fluxes are in the reflux.  Tests/CTestList.cmake runs
# it on several meshes, with map factors, and with a sheared flow and turbulent diffusion; and,
# with a uniform wind over a hill, checks that the composite mass ("MASS SL/ML") is conserved.
#
erf.prob_name = "Scalar Advection/Diffusion"

max_step = 50

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 1  1  1
amr.n_cell           = 32 32 32

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt       = 5e-5
erf.fixed_fast_dt  = 2.5e-5

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1
erf.v              = 1
amr.v              = 1

# REFINEMENT -- given as a vector, so that a test can refine anisotropically
amr.max_level = 1
amr.ref_ratio_vect = 2 2 2
erf.refinement_indicators = box1
erf.box1.max_level = 1
erf.box1.in_box_lo = 0.25 0.25 0.0
erf.box1.in_box_hi = 0.75 0.75 0.5
erf.coupling_type  = "TwoWay"
erf.cf_width       = 0
erf.cf_set_width   = 0

# CHECKPOINT FILES
erf.check_file      = chk
erf.check_int       = -1

# PLOTFILES
erf.plot_file_1     = plt
erf.plot_int_1      = -1
erf.plot_vars_1     = density rhoadv_0

# SOLVER CHOICE
erf.use_gravity = false

erf.les_type        = "None"
erf.molec_diff_type = "Constant"
erf.alpha_T         = 0.0
erf.alpha_C         = 1.0

erf.dycore_horiz_adv_type  = "Centered_2nd"
erf.dycore_vert_adv_type   = "Centered_2nd"
erf.dryscal_horiz_adv_type = "Centered_2nd"
erf.dryscal_vert_adv_type  = "Centered_2nd"

erf.init_type = "uniform"

# PROBLEM PARAMETERS -- a blob (prob_type 10) at rest
prob.rho_0 = 1.0
prob.A_0 = 1.0
prob.U_0 = 0.0
prob.V_0 = 0.0
prob.prob_type = 10
