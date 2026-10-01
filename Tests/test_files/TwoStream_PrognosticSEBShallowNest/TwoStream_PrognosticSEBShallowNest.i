# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# A SHALLOW nest -- a fine level no box of which spans the domain in z. Such a level
# cannot sweep: TwoStreamRadiation::advance returns at its top and advance_radiation
# interpolates its radiation fields from the parent instead. Its surface energy balance
# never runs either, so its t_sfc/q_sfc stay frozen at whatever fill_seb_from_coarse
# wrote when the level was built.
#
# The danger is the average-down. Averaging that frozen field onto level 0 would pin the
# coarse surface under the patch at its level-creation value, so the longwave boundary
# condition there stops responding to the surface energy balance -- a regression against
# the old level-0-only behaviour. post_timestep therefore skips the transfer for a level
# that does not sweep (ERF::rad_level_needs_interpolation).
#
# This is the complement of TwoStream_PrognosticSEBMultiLevel, whose
# amr.refine_whole_domain_dir = 2 makes its fine boxes span z. Neither case covers the
# other: that one checks the transfer happens, this one checks it does not.
#
# The grid settings come from TwoStream_NestedPatch, where the shallow patch is already
# established (k = 0..7 of the 32-cell domain). Do not add terrain or change value_less
# without re-checking the patch is still shallow -- if level 1 ends up spanning z the
# guard never engages and this case proves nothing.
#
# The checker asserts every cell of level 0 moved away from erf.rad_t_sfc; with the guard
# removed, the cells under the nest sit exactly at it.

erf.prob_name = "ABL"

max_step = 6
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

# No amr.refine_whole_domain_dir here: tagging on theta (which increases with
# height) tags only the lower part of the column, so level 1 comes out as a
# shallow patch -- k = 0..7 of the 32-cell domain. That is a nested patch: it
# carries no complete column, so the sweep cannot run on it and its heating
# rates are interpolated from level 0 instead.
erf.refinement_indicators = lowth
erf.lowth.max_level = 1
erf.lowth.field_name = theta
erf.lowth.value_less = 301.0

erf.check_file = chk
erf.check_int = -1

erf.plot_file_1 = plt
erf.plot_int_1 = 2
erf.plot_vars_1 = density theta qsrc_sw qsrc_lw

erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 2
erf.plot2d_vars_1 = seb_t_sfc seb_q_sfc

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
# Drive the surface from the sweep's own fluxes, or T_s never leaves erf.rad_t_sfc and
# the check passes on a constant.
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_restore_timescale_s = 3600.0
