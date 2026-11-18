# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Line coordinates of a level made by a regrid after the start: amr.max_level = 1, but the
# refinement criterion (theta > 0 everywhere) only starts at t = 0.05, so level 1 does not
# exist at the start and covers the domain from the regrid after step 1. The text output of
# a line given by real coordinates must then write plt_line_lev1.zcc, the heights of the
# column on level 1 (check_line_zcc.py).
erf.prob_name = "ABL"
max_step = 3

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 256 256 256
amr.n_cell           =  16  16  16
amr.max_grid_size    =   8
geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 0.1

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.v            = 1

# REFINEMENT: level 1 is made by the regrid after step 1, over the whole domain
amr.max_level = 1
amr.ref_ratio = 2
erf.regrid_int = 1
erf.refinement_indicators = grow
erf.grow.max_level = 1
erf.grow.field_name = theta
erf.grow.value_greater = 0.0
erf.grow.start_time = 0.05

# LINE SAMPLING: one vertical line through the column at (56, 72), written as text
erf.do_line_sampling = true
erf.line_sampling_interval = 1
erf.line_sampling_text_output = true
erf.sample_line_lo_real = 56.0 72.0 1.0
erf.sample_line_hi_real = 56.0 72.0 255.0
erf.sample_line_dir = 2
erf.sample_line_name = col
erf.line_sampling_vars = x_velocity

# CHECKPOINT AND PLOTFILES
erf.check_int  = -1
erf.plot_int_1 = -1

# SOLVER CHOICE
erf.use_gravity     = false
erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.1
# The fine grids are split in z, which implicit vertical diffusion does not allow
erf.vert_implicit_fac = 0 0 0

# PROBLEM PARAMETERS
erf.init_type = "uniform"
prob.rho_0 = 1.0
prob.A_0   = 1.0
prob.U_0   = 10.0
prob.V_0   = 0.0
prob.W_0   = 0.0
prob.T_0   = 300.0
prob.U_0_Pert_Mag = 0.08
prob.V_0_Pert_Mag = 0.08
prob.W_0_Pert_Mag = 0.0
