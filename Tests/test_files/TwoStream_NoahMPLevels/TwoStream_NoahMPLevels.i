# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# Two-stream radiation feeding Noah-MP on two levels. The case is
# Exec/RegTests/NoahMP_Ideal/inputs_noahmp_twostream (a 4 x 4 grassland patch at
# 250 m, 18:00 UTC at 40 N, 100 W, one 3600 s land step on the first ERF step) with a
# level 1 over the middle 2 x 2 cells, refined by 2 in x and y. Needs ERF_ENABLE_NOAHMP.
#
# The runner (Tests/RunTwoStreamNoahMPLevels.cmake) runs it in four ways:
#
#   own     namelist.erf names ERF_SETUP_FILE_02 = wrfinput_d02, a barren nested land
#           file, so level 1 runs Noah-MP itself, on the forcing of its own sweep
#           (amr.refine_whole_domain_dir = 2: the level spans z and sweeps its columns);
#   nested  the same land, but the level stops below the domain top, so it does not
#           sweep: its radiation, Noah-MP's forcing included, comes from level 0;
#   interp  no ERF_SETUP_FILE_02: level 1 takes its land state from level 0;
#   regrid  the refined region moves at t = 0.5 s; a level that runs Noah-MP on its own
#           land file cannot be rebuilt, so the run must stop.
#
# check_two_stream_noahmp_levels.py asserts what each must show; see its docstring.
erf.prob_name = "ABL"

max_step = 2
stop_time = 100.0
amrex.fpe_trap_invalid = 0

geometry.prob_extent = 1000 1000 1000
amr.n_cell           = 4 4 32
geometry.is_periodic = 1 1 0

amr.max_level = 1
amr.ref_ratio_vect = 2 2 1
amr.n_error_buf = 0
amr.blocking_factor = 2
amr.max_grid_size = 64
# Every level spans z and sweeps its own columns (the nested leg turns this off).
amr.refine_whole_domain_dir = 2
# Level 1 over the middle 2 x 2 coarse cells: parent cell (2, 2), 1-based, is the first,
# as I_PARENT_START and J_PARENT_START in wrfinput_d02 say.
erf.refinement_indicators = patch
erf.patch.max_level = 1
erf.patch.in_box_lo = 250.0 250.0
erf.patch.in_box_hi = 750.0 750.0
# No subcycling: each level sweeps once per step, so the step-1 row of the radiation
# CSV is the sweep the land step integrated on.
erf.dt_ref_ratio = 1

zlo.type = "surface_layer"
zhi.type = "SlipWall"
erf.most.z0   = 0.1
erf.most.zref = 20.0

erf.fixed_dt = 1.0
erf.v = 1

erf.init_type = "input_sounding"
erf.sounding_type = Ideal
erf.input_sounding_file = "input_sounding"
erf.use_gravity = true
erf.les_type = "Smagorinsky"
erf.Cs = 0.16
erf.pbl_type = "None"
erf.molec_diff_type = "None"

start_datetime = "2024-08-05 18:00:00"
erf.rad_cons_lat = 40.0
erf.rad_cons_lon = -100.0

erf.radiation_model = "TwoStream"
erf.rad_t_sfc = 300.0
erf.radiation.tau_per_layer    = 0.003125
erf.radiation.tau_lw_per_layer = 0.05
erf.radiation.diag_enable = true
erf.radiation.diag_csv_enable = true
erf.radiation.diag_stdout_enable = false
erf.radiation.diag_callsite_mode = pre_only
erf.radiation.diag_file = "radiation_diag.csv"

erf.land_surface_model = "NOAHMP"

erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1 = -1
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 2
erf.plot2d_vars_1 = t_sfc sav sag albedo sw_flux_dn lw_flux_dn cos_zenith_angle
