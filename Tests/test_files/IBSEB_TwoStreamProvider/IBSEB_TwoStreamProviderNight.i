# IBSEB_TwoStreamProvider.i at night under the two-stream calendar sun: 06:00 UTC on
# 5 August at 40 N, 100 W, about 23:20 local solar time, with the faces on the same sun
# (erf.ibseb.sun_mode = two_stream), an absorbing sky, floating-point traps on, and a 70 m
# tower (20 m x 20 m at x, y = 40-60 m) beside the 40 m cube. The faces must get no
# shortwave and the sky's longwave, less of it on the higher roof, and nothing may divide
# by the zero cosine of the zenith. Run by Tests/RunIBSEBTwoStreamProvider.cmake.
stop_time = 100000.0
max_step  = 1
amrex.fpe_trap_invalid = 1
amrex.fpe_trap_zero = 1
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
erf.buildings_file_name = cube_and_tower_32x32.txt
erf.immersed_forcing_substep = true
eb2.small_volfrac = 0.005
erf.if_use_most = true
erf.if_z0 = 0.01
erf.use_coriolis = false

erf.ibseb.enable = true
erf.ibseb.prognostic = false
erf.ibseb.debug = false
erf.ibseb.csv_int = 1
erf.ibseb.csv_file = "ibseb_buildings.csv"
erf.ibseb.dump_faces_file = "faces"
erf.ibseb.dump_faces_tag_step = true
erf.ibseb.T_skin_init = 300.0
erf.ibseb.T_interior = 293.0
erf.ibseb.n_slab_layers = 4
erf.ibseb.albedo = 0.3
erf.ibseb.emissivity = 0.9
erf.ibseb.radiation = "two_stream"
erf.ibseb.sun_mode = "two_stream"

erf.radiation_model = "TwoStream"
start_datetime = "2024-08-05 06:00:00"
erf.rad_cons_lat = 40.0
erf.rad_cons_lon = -100.0
erf.rad_t_sfc = 300.0
erf.radiation.tau_per_layer = 0.02
erf.radiation.tau_lw_per_layer = 0.3
erf.radiation.surface_albedo_sw = 0.2
erf.radiation.surface_emissivity_lw = 0.95
