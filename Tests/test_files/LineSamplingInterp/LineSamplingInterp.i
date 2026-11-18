# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Lines interpolated to requested heights (erf.sample_line_interp) on a terrain-fitted mesh.
# Over a 100 m Witch of Agnesi hill, on 12 levels from 20 m stretched by 1.2, a scalar that
# varies in x and z (prob.prob_type = 12) is sampled on five lines:
#   tower      over the hill top, 15 heights from 0 to 280 m
#   far        far from the hill, 15 heights from 0 to 280 m
#   high       over the hill top, 5 heights from 600 to 790 m, the highest above the top cell
#              centre; above the ground (694 to 884 m) one is bracketed, two lie between the top
#              cell centre and the 791.6 m top, and two above the top, which are missing
#   one        over the hill top at the single height 150 m
#   across     in x at y = 40 m across the hill at z = 50 m, part of it below the ground
# The CTest entry runs it as given (heights above sea level) and with erf.sample_line_z_agl = 1
# (heights above the ground under each column), and check_line_interp.py interpolates the
# plotfile's rhoadv_0 to the same heights and compares.  LineSamplingInterp_TwoLevel runs the
# same deck with a full-depth refined level over the hill.

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

# LINE SAMPLING: x, y = 40 is cell (16,16), over the hill top
erf.do_line_sampling          = true
erf.line_sampling_interval    = 1
erf.line_sampling_text_output = true
erf.line_sampling_vars        = conserved
erf.sample_line_interp        = true
erf.sample_line_dir     = 2 2 2 2 0
erf.sample_line_npts    = 15 15 5 1 0
erf.sample_line_name    = tower far high one across
erf.sample_line_lo_real =   40.   40.   0.   -1240. -1240.   0.   40. 40. 600.   40. 40. 150.   -1240. 40. 50.
erf.sample_line_hi_real =   40.   40. 280.   -1240. -1240. 280.   40. 40. 790.   40. 40. 150.    1240. 40. 50.

# CHECKPOINT AND PLOTFILES
erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 1
erf.plot_vars_1 = rhoadv_0 z_phys

# SOLVER CHOICE
erf.dycore_horiz_adv_type  = "Upwind_3rd"
erf.dycore_vert_adv_type   = "Upwind_3rd"
erf.dryscal_horiz_adv_type = "Upwind_3rd"
erf.dryscal_vert_adv_type  = "Upwind_3rd"

erf.use_gravity     = false
erf.les_type        = "None"
erf.molec_diff_type = "None"

# PROBLEM PARAMETERS: a scalar A_0 (1 + cos(pi r / r0)) / 4 in x and z, r0 = 0.25 of the domain width
erf.init_type = "uniform"
prob.rho_0     = 1.0
prob.T_0       = 300.0
prob.A_0       = 1.0
prob.U_0       = 10.0
prob.V_0       = 3.0
prob.prob_type = 12
prob.rad_0     = 0.25
