# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# qsrc_sw and qsrc_lw must be advertised as available plotfile variables
# whenever qheating_rates is allocated, which includes the TwoStream
# radiation path (not only RRTMGP). This zero-step run checks that both
# names survive plot-variable selection and appear in the plotfile header.
erf.prob_name = "ABL"

max_step = 0
stop_time = 1.0
amrex.fpe_trap_invalid = 0

fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_extent = 1024 1024 1024
amr.n_cell           = 4 4 16
amrex.max_grid_size_z = 128
geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"
zhi.theta_grad = 0.003

erf.fixed_dt = 0.5
erf.sum_interval = 1
erf.v = 1
amr.v = 1
amr.max_level = 0

erf.check_file = chk
erf.check_int = -1

erf.plot_file_1 = plt
erf.plot_int_1 = 1
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

erf.radiation_type = "TwoStream"
erf.radiation.sw_enabled = true
erf.radiation.lw_enabled = true
erf.radiation.v = 0
erf.radiation.diag_csv_enable = false
erf.radiation.diag_stdout_enable = false
