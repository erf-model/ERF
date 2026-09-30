# Two levels initialised from an input sounding (issue #4143).
#
# Every level samples the sounding, and the large-scale forcing profiles, at its own cell
# centres. The soundings have an entry every 5 m, on every level-0 (10, 30, ... m) and level-1
# (5, 15, ... m) cell centre, with values exact in binary: theta = 300 + z/40 and U = z/64 in
# input_sounding_sfc300/303, theta = 300 and U = z/64 in input_sounding_wind. Four tests use
# this deck:
#
# - InputSounding_FineLevelInit: the two sfc files differ only in their surface line (303 K
#   against 300 K). A fine cell samples the file at its own height, above that line, so both
#   runs start from the same state. Resampling a level-0 profile instead blended the surface
#   line into the fine cells below the first level-0 centre (the lowest fine layer then
#   starts at 301.625 K with a 303 K surface line).
#
# - InputSounding_FineLevelNudging (theta) and InputSounding_FineLevelWindNudging (u and v,
#   through the momentum sources): one step with and without nudging towards the sounding the
#   run started from, so nudging adds exactly zero on every level. The level-0 profile indexed
#   by a fine k is the value from twice the height, which nudged the fine level.
#
# - InputSounding_FineLevelLSF: one step with and without large-scale forcing whose
#   tendencies and subsidence are zero and whose wind equals the sounding's
#   (lsf_zero_tendency), so relaxing towards it adds exactly zero on every level.

erf.prob_name = "ABL"
max_step = 0
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo     = 0.0 0.0 0.0
geometry.prob_hi     = 80.0 80.0 320.0
amr.n_cell           = 4 4 16
amr.max_grid_size    = 256
amr.blocking_factor  = 2
geometry.is_periodic = 1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"

amr.max_level        = 1
amr.ref_ratio_vect   = 2 2 2
erf.dt_ref_ratio     = 1
erf.coupling_type    = "TwoWay"
erf.refinement_indicators = box1
erf.box1.max_level = 1
erf.box1.in_box_lo_indices_crse = 0 0 0
erf.box1.in_box_hi_indices_crse = 3 3 7

erf.fixed_dt = 0.01
erf.init_type           = "input_sounding"
erf.input_sounding_file = "../input_sounding_sfc300"
erf.sounding_type       = "ConstantDensity"
erf.use_gravity = false
erf.molec_diff_type = "None"
erf.les_type        = "None"

erf.v = 1
amr.v = 0
erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 1
erf.plot_vars_1 = density theta x_velocity z_velocity
