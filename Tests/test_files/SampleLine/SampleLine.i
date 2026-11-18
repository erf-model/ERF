# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# erf.sample_line is a LineSampler of full-depth columns on level 0.  This deck samples two
# columns both ways, with erf.sample_line (one log per column) and with erf.do_line_sampling
# (text output, one file per variable) with the variables erf.sample_line writes;
# check_sample_line.py requires the two to agree at every sample.  The boxes are 8 x 8, so
# the columns are on different boxes and, on 2 ranks, on different ranks.

erf.prob_name = "ABL"

max_step = 3

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 256 256 256
amr.n_cell           =  16  16  16
amr.max_grid_size    =   8
amr.max_level        =   0
geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 0.1

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.v            = 1

# Two columns: the k of each cell only selects the column
erf.sample_line_log = col_a col_b
erf.sample_line     = 3 4 0  12 10 7

erf.do_line_sampling           = true
erf.line_sampling_interval     = 1
erf.line_sampling_text_output  = true
erf.sample_line_lo   = 3 4 0   12 10 0
erf.sample_line_hi   = 3 4 15  12 10 15
erf.sample_line_dir  = 2 2
erf.sample_line_name = dls_a dls_b
erf.line_sampling_vars = conserved x_velocity y_velocity z_velocity sgs_tau11 sgs_tau12 sgs_tau13 sgs_tau22 sgs_tau23 sgs_tau33

# CHECKPOINT AND PLOTFILES
erf.check_int  = -1
erf.plot_int_1 = -1

# SOLVER CHOICE
erf.use_gravity     = false
erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.1

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
