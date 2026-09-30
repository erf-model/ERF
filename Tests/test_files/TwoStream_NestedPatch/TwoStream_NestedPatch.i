# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# A nested patch -- a fine level that stops below the domain top -- is a
# supported configuration, not an error. ERF interpolates its heating rates and
# fluxes from the parent, the same route RRTMGP takes (is_nested_patch).
#
# The runner checks levels 0 and 1. On level 1 the checker detects the nested
# patch from the data and asserts what remains meaningful there: finite,
# non-negative, and not identically zero -- which is exactly what fails if the
# interpolation never happens and qheating_rates keeps the zeros it was
# allocated with. The column-structure assertions are skipped there because a
# patch that does not contain the top of the atmosphere has no such layer.

erf.prob_name = "ABL"

max_step = 2
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
erf.radiation.diag_csv_enable = false
erf.radiation.diag_stdout_enable = false
