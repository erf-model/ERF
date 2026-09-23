# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# TwoStream radiation on a refined hierarchy. Every level runs its own column
# sweep, so the vertical structure the single-level case checks must hold on
# the fine level too; the runner invokes TwoStreamRadiationCheck with
# CHECK_LEVELS 0 1.
#
# What this catches that the single-level cases cannot:
#  - a fine level that is never swept. qheating_rates[lev] is allocated and
#    setVal(zero) for every level, and the RhoTheta source applies it with no
#    level gate, so a missed fine level is silently zero heating -- which fails
#    the checker's "qsrc_sw is zero everywhere" assertion at level 1.
#  - the refinement patch is tagged on a field, not an explicit erf.boxN, so it
#    exercises the path ERF_RefineBox.cpp does not cover. Without
#    amr.refine_whole_domain_dir the tagged patch stops short of the domain top
#    and TwoStreamRadiation::define_level aborts; with it the patch spans z.
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

# Tag on theta, which increases with height (zhi.theta_grad above): the tagged
# region is the lower part of the column and would give a patch that stops
# short of the domain top. refine_whole_domain_dir = 2 makes AMReX cluster in
# the horizontal only and emit boxes that span z, which is what the column
# sweep requires (and what the model's abort message recommends).
amr.refine_whole_domain_dir = 2
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
