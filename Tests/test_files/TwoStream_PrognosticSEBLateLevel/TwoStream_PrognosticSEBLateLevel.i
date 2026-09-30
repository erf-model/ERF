# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# A fine level created MID-RUN must start from the surface its parent has
# reached. Every level evolves its own force-restore surface temperature, and a
# level built by a regrid gets it from ERF::fill_seb_from_coarse, which
# interpolates the parent's t_sfc/q_sfc. Without that the new level starts from
# the erf.rad_t_sfc scalar -- and the average-down then drags the coarse surface
# under it back there as well.
#
# TwoStream_PrognosticSEBMultiLevel cannot see this: its fine level exists from
# t = 0, when both levels are still uniform at erf.rad_t_sfc, so interpolating
# from the parent and filling with the scalar give the same field. Here the
# refinement criterion only switches on at erf.lowth.start_time, so level 0 runs
# alone for 10 steps and its surface drifts ~0.08 K (with x-structure from the
# hill) before the regrid at the start of step 11 builds level 1.
#
# check_two_stream_seb_parity.py --created-from asserts that level 1 at step 11,
# one step after its creation, averages over each coarse cell to the parent's
# surface extrapolated one step past step 10 -- to a quarter of one step's
# change. A level started from the scalar misses by the full 10-step drift.
#
# Everything else is TwoStream_PrognosticSEBMultiLevel's deck; keep the two in
# step.

erf.prob_name = "ABL"

max_step = 11
stop_time = 10.0
amrex.fpe_trap_invalid = 0

# Deliberately left on AMReX's default MFIter tile size, which splits the
# domain in z. The column sweep must be independent of that tiling; an
# earlier version restarted the sweep at the bottom of every z tile and
# this case is what catches that.

geometry.prob_extent = 1024 1024 1024
amr.n_cell           = 4 4 32
amr.max_grid_size_z = 128
geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"
zhi.theta_grad = 0.003

erf.fixed_dt = 0.5
erf.sum_interval = 1
erf.v = 1
amr.v = 1
amr.max_level = 1
amr.ref_ratio_vect = 2 2 1
amr.n_error_buf = 0

# Tag on theta, which increases with height (zhi.theta_grad above): the tagged
# region is the lower part of the column and would give a patch that stops
# short of the domain top. refine_whole_domain_dir = 2 makes AMReX cluster in
# the horizontal only and emit boxes that span z, which is what the column
# sweep requires (and what the model's abort message recommends).
#
# start_time holds the tagging off until t = 4.75 s, so with regrid_int = 1 the
# first regrid that sees a tag is the one at the start of step 11 (t = 5.0 s).
# The half-step margin keeps that independent of how t = 5.0 rounds.
amr.refine_whole_domain_dir = 2
erf.refinement_indicators = lowth
erf.lowth.max_level = 1
erf.lowth.field_name = theta
erf.lowth.value_less = 1.0e4   # above every theta: tags the whole domain
erf.lowth.start_time = 4.75
erf.regrid_int = 1

erf.check_file = chk
erf.check_int = -1

erf.plot_file_1 = plt
erf.plot_int_1 = 11
erf.plot_vars_1 = density theta qsrc_sw qsrc_lw

# The surface state itself, which is what the check reads: every step, so it
# has plt2d00009 and plt2d00010 (level 0 alone) and plt2d00011 (level 1 new).
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 1
erf.plot2d_vars_1 = seb_t_sfc seb_q_sfc

# Terrain: a Witch of Agnesi hill on a fitted mesh, so every column has its own
# layer thicknesses and the fine level resolves the hill the coarse one cannot.
# Without this the problem is horizontally uniform, both levels compute identical
# surface fluxes and evolve identically -- and the test would pass whether or not
# the levels are actually kept in step, which is worth nothing.
erf.terrain_type         = StaticFittedMesh
erf.terrain_smoothing    = 0
prob.custom_terrain_type = "WoA"
prob.dir                 = 0
prob.hmax                = 100.0
prob.L                   = 300.0

erf.use_gravity = true
erf.molec_diff_type = "None"
erf.les_type = "None"
erf.pbl_type = "None"
erf.theta_ref = 300.0

erf.init_type = "input_sounding"
erf.sounding_type = Ideal
erf.input_sounding_file = "input_sounding"

erf.use_coriolis = false
erf.abl_driver_type = "None"

# RADIATION - TwoStream, SW + LW, clear sky, fixed sun
erf.radiation_model = "TwoStream"
erf.radiation.sw_enabled = true
erf.radiation.lw_enabled = true
erf.radiation.tau_per_layer = 0.00625
erf.radiation.tau_lw_per_layer = 1.0
erf.fixed_solar_zenith_angle = 0.5    # cos(60 deg): the cosine, as RRTMGP takes it
erf.fixed_total_solar_irradiance = 1361.0
erf.rad_t_sfc = 300.0    # surface temperature [K] where no LSM or surface layer supplies one (shared with RRTMGP)
erf.radiation.v = 0
erf.radiation.diag_csv_enable = true
erf.radiation.diag_stdout_enable = true

# The prognostic surface energy balance: the feature under test.
erf.radiation.seb_enable = true
erf.radiation.seb_prognostic_enable = true
erf.radiation.diag_enable = true

# The surface has to actually move, or the test passes on a constant and proves
# nothing: with the defaults every flux in the balance is zero and T_s sits at
# erf.rad_t_sfc forever. Drive it from the sweep's own surface fluxes and shorten
# the restore timescale so the drift is resolvable in a few steps.
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_restore_timescale_s = 3600.0
