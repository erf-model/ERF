# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Vertical lines given by real coordinates on a terrain-fitted mesh: their cells come from the
# heights of each column, and the height of each sample is written with it (the variable z).
# Over a 100 m Witch of Agnesi hill, on 12 levels from 20 m stretched by 1.2, four lines are
# sampled: three vertical ones, over the hill top and far from it from z = 0 to 300 m and over
# the hill top from 150 to 420 m, and one horizontal one across the hill at z = 150 m.  The
# CTest entry runs it as given (heights above sea level) and with erf.sample_line_z_agl = 1
# (heights above the ground), and check_line_heights.py compares each line with z_phys in the
# plotfile.

erf.prob_name = "Scalar Advection/Diffusion"
erf.anelastic = 1

max_step  = 1
stop_time = 1.0e6

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -1280. -1280.    0.
geometry.prob_hi     =  1280.  1280.  791.6
amr.n_cell           =    32     32   12
geometry.is_periodic = 1 1 0
amr.max_level        = 0

erf.terrain_type          = StaticFittedMesh
erf.terrain_smoothing     = 0
erf.initial_dz            = 20.0
erf.grid_stretching_ratio = 1.2
prob.custom_terrain_type  = "WoA"
prob.dir                  = 2
prob.hmax                 = 100.0
prob.L                    = 300.0
erf.terrain_poisson_solver = mlmg

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 2.0

# DIAGNOSTICS & VERBOSITY
erf.v = 1
amr.v = 0

# LINE SAMPLING: cells (16,16) and (0,0) are at x, y = 40 and -1240
erf.do_line_sampling          = true
erf.line_sampling_interval    = 1
erf.line_sampling_text_output = true
erf.line_sampling_vars        = theta
# and the fourth line runs in x at y = 40 over the hill, on the level of cells that holds
# z = 150 m at its first point (cell 0,16)
erf.sample_line_dir     = 2 2 2 0
erf.sample_line_name    = top far top_high across
erf.sample_line_lo_real =   40.   40.   0.   -1240. -1240.   0.   40. 40. 150.   -1240. 40. 150.
erf.sample_line_hi_real =   40.   40. 300.   -1240. -1240. 300.   40. 40. 420.    1240. 40. 150.

# CHECKPOINT AND PLOTFILES
erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 1
erf.plot_vars_1 = theta z_phys

# SOLVER CHOICE
erf.dycore_horiz_adv_type  = "Upwind_3rd"
erf.dycore_vert_adv_type   = "Upwind_3rd"
erf.dryscal_horiz_adv_type = "Upwind_3rd"
erf.dryscal_vert_adv_type  = "Upwind_3rd"

erf.use_gravity     = false
erf.les_type        = "None"
erf.molec_diff_type = "None"

# PROBLEM PARAMETERS
erf.init_type = "uniform"
prob.rho_0 = 1.0
prob.T_0   = 300.0
prob.A_0   = 1.0
prob.U_0   = 10.0
prob.V_0   = 3.0
prob.rad_0 = 100.0
