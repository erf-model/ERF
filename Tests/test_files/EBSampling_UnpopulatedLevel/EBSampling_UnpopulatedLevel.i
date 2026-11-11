# Line and plane sampling with EB terrain (#4977) while a fine level allowed by
# amr.max_level does not exist: the refinement criterion never triggers, so level 1
# is never made and its EB factory is never built.  The column through the top of the
# Witch of Agnesi hill and a plane through its flanks must write covered cells as
# erf.sampling_covered_value (Tests/test_files/EBSampling_UnpopulatedLevel/check_eb_sampling.py)
erf.prob_name = "Flow over Witch of Agnesi hill"
max_step = 2
erf.init_type = ConstantDensity

amrex.fpe_trap_invalid=1
amrex.fpe_trap_zero=1
amrex.fpe_trap_overflow=1

fabarray.mfiter_tile_size = 1024 1024 1024

eb2.geometry = terrain
eb2.small_volfrac = 1.e-4

geometry.prob_lo     = 0.   0.   0.
geometry.prob_hi     = 595.0  20.0  600.
amr.n_cell           = 60  4   86
geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.timestepping_type = RK3_NoSubstepping
erf.cfl              = 0.5

erf.sum_interval   = 1
erf.v              = 1
amr.v              = 1

# A fine level is allowed but never made
amr.max_level = 1
amr.ref_ratio = 2
erf.regrid_int = 1
erf.refinement_indicators = never
erf.never.max_level = 1
erf.never.field_name = theta
erf.never.value_greater = 1.0e10

erf.check_int      = -1
erf.plot_int_1     = -1

erf.use_gravity = true
erf.use_coriolis = false
erf.les_type = "None"

erf.terrain_type = "EB"
erf.terrain_smoothing = 0
erf.eb_boundary_type = "NoSlipWall"

erf.molec_diff_type   = "Constant"
erf.dynamic_viscosity = 60.0
erf.alpha_T           = 0.0
erf.theta_ref = 300.0
erf.eb_diff_constraint_y = true

erf.abl_driver_type = "PressureGradient"
erf.abl_pressure_grad = -0.005 0. 0.

erf.project_initial_velocity = 0

prob.custom_terrain_type = "WoA"
prob.dir   = 0
prob.T_0   = 300.0
prob.U_0   = 0.0
prob.V_0   = 0.0
prob.rho_0 = 1.16
prob.hmax = 100.0
prob.L = 100.0

# The full-depth column through the cell at the top of the hill (i = 30), written to a log
erf.sample_line_real = 297.5 10.0 0.0
erf.sample_line_log  = col_top

# The same column with erf.do_line_sampling, written as text (which also writes the
# heights of each level, plt_line_lev*.zcc)
erf.do_line_sampling = true
erf.line_sampling_interval = 1
erf.line_sampling_text_output = true
erf.sample_line_lo_real = 297.5 10.0 1.0
erf.sample_line_hi_real = 297.5 10.0 599.0
erf.sample_line_dir = 2
erf.sample_line_name = line_top
erf.line_sampling_vars = volfrac density theta

# A plane through both flanks of the hill
erf.do_plane_sampling = true
erf.plane_sampling_interval = 1
erf.sample_plane_lo = 0.0 0.0 30.0
erf.sample_plane_hi = 595.0 20.0 30.0
erf.sample_plane_dir = 2
erf.sample_plane_name = plane_z30
erf.plane_sampling_vars = density theta volfrac
