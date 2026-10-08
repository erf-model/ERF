# The building faces on the two-stream columns (erf.ibseb.radiation = two_stream): the
# 40 m cube of Tests/test_files/IBSEB_Cube under the two-stream radiation with a fixed
# sun (cos z 0.766044443118978, 40 degrees, from the south-south-west for the faces) and
# a transparent sky (no shortwave or longwave optical depth), the skin held at its
# initial temperature so each term can be compared. IBSEB_TwoStreamProvider
# (Tests/RunIBSEBTwoStreamProvider.cmake) runs it as given, again with the prescribed
# provider set to the same sky (direct-normal irradiance 1000 W/m2, no diffuse light,
# no sky longwave, the ground's albedo, emissivity and temperature), and again with an
# absorbing sky, and checks the faces with check_ibseb_two_stream_provider.py.
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
erf.ibseb.sun_mode = "fixed"
erf.ibseb.sun_zenith_deg = 40.0
erf.ibseb.sun_azimuth_deg = 200.0

erf.radiation_model = "TwoStream"
erf.fixed_solar_zenith_angle = 0.766044443118978
erf.fixed_total_solar_irradiance = 1000.0
erf.rad_t_sfc = 300.0
erf.radiation.tau_per_layer = 0.0
erf.radiation.tau_lw_per_layer = 0.0
erf.radiation.surface_albedo_sw = 0.2
erf.radiation.surface_emissivity_lw = 0.95
