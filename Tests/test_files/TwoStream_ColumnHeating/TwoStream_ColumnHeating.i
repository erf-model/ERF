# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# The TwoStream radiation sweep must treat k = 0 as the surface layer and
# the highest k as the top of the atmosphere, and the LW heating sign must
# give cooling to space. This short SW + LW run writes the per-level heating
# rates (qsrc_sw, qsrc_lw) to plt00002; TwoStreamRadiationCheck verifies the
# vertical structure: SW heating strongest at the top layer, LW cooling
# strongest at the top layer, net LW cooling of the column.
erf.prob_name = "ABL"

max_step = 2
stop_time = 10.0
amrex.fpe_trap_invalid = 0

fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_extent = 1024 1024 1024
amr.n_cell           = 4 4 32
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
erf.radiation_type = "TwoStream"
erf.radiation.sw_enabled = true
erf.radiation.lw_enabled = true
erf.radiation.tau_per_layer = 0.00625
erf.radiation.tau_lw_per_layer = 1.0
erf.radiation.solar_zenith = 60.0
erf.radiation.S0 = 1361.0
erf.radiation.surface_temp_k = 300.0
erf.radiation.v = 0
erf.radiation.diag_csv_enable = false
erf.radiation.diag_stdout_enable = false
