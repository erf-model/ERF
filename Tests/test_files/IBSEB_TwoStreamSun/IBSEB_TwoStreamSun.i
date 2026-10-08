# The face balance under the two-stream radiation: the 40 m cube of
# Tests/test_files/IBSEB_Cube with erf.radiation_model = TwoStream following the
# calendar, and the faces on the two-stream sun (erf.ibseb.sun_mode = two_stream),
# which ERF::ibseb_check_sun_matches_two_stream() requires. IBSEB_TwoStreamSunRun
# runs it and compares the faces' sun with the sweep's diagnostics; the abort
# tests break it on the command line (erf.ibseb.sun_mode = solar, a fixed
# two-stream sun, no two-stream radiation).
stop_time = 100000.0
max_step  = 1
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_extent =  320     320    160
amr.n_cell           =  32      32     16
amr.max_grid_size    =  16
geometry.is_periodic = 1 1 0
zlo.type = "NoSlipWall"
zhi.type = "SlipWall"
erf.fixed_dt           = 0.5
erf.substepping_cfl    = 0.5
erf.sum_interval   = -1
erf.v              = 0
amr.v              = 0
amr.max_level      = 0
erf.check_file      = chk
erf.check_int       = -1
erf.plot_file_1     = plt
erf.plot_int_1      = -1
erf.plot_vars_1     = theta terrain_IB_mask ibseb_nfaces ibseb_tskin ibseb_sw_abs ibseb_lw_net ibseb_H ibseb_G
erf.use_gravity = true
erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.17
erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.buildings_type = ImmersedForcing
erf.buildings_file_name = cube_40m_10m_32x32.txt
erf.immersed_forcing_substep = true
eb2.small_volfrac = 0.005
erf.if_use_most = true
erf.if_z0 = 0.01
erf.use_coriolis = false

erf.ibseb.enable = true
erf.ibseb.prognostic = true
erf.ibseb.debug = false
erf.ibseb.csv_int = -1
erf.ibseb.T_skin_init = 300.0
erf.ibseb.T_interior = 293.0
erf.ibseb.n_slab_layers = 8
erf.ibseb.k_therm = 1.0
erf.ibseb.rho_cp = 1.6e6
erf.ibseb.thickness = 0.2
erf.ibseb.albedo = 0.3
erf.ibseb.emissivity = 0.9
erf.ibseb.z0_wall = 0.01
erf.ibseb.z0h_wall = 0.001
erf.ibseb.lw_mode = "gray"
erf.ibseb.sky_emissivity = 0.83
erf.ibseb.T_ground = 300.0
erf.ibseb.newton_tol_K = 1.0e-3
erf.ibseb.newton_max_iter = 20

erf.radiation_model = "TwoStream"
erf.rad_t_sfc = 300.0
start_datetime = "2024-08-05 15:00:00"
erf.rad_cons_lat = 40.0
erf.rad_cons_lon = -100.0
erf.ibseb.sun_mode = "two_stream"
erf.ibseb.sw_transmission = 0.8
erf.radiation.diag_enable = true
erf.radiation.diag_csv_enable = true
erf.radiation.diag_stdout_enable = false
erf.radiation.diag_callsite_mode = pre_only
erf.radiation.diag_file = "radiation_diag.csv"
erf.ibseb.csv_file = "ibseb_buildings.csv"
