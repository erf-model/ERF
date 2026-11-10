# The two-stream columns and their ground balance seeing the buildings
# (erf.ibseb.two_stream_buildings, on by default): the deck of
# Tests/test_files/IBSEB_TwoStreamProvider (the 40 m cube under a transparent sky, the sun
# fixed 40 degrees from the zenith in the south-south-west, the faces on the two-stream
# columns) with the force-restore ground balance on the sweep's shortwave, and the ground
# and heating fields in the plotfiles. IBSEB_TwoStreamBuildings
# (Tests/RunIBSEBTwoStreamBuildings.cmake) runs it as given, restarted, under an absorbing
# and a scattering sky, with the input off, at night and on two levels, and checks the
# fields with check_ibseb_two_stream_buildings.py.
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
erf.plot_int_1      = 1
erf.plot_vars_1     = theta terrain_IB_mask qsrc_sw qsrc_lw
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
# The roof (z = 40 m) lies on a plane of the grid refined by 2 that the face-centred
# EB factories are built on; in single precision that gives one face with more than
# two cuts (erf-model/ERF#5006). In double precision there are none, so this changes nothing.
eb2.cover_multiple_cuts = 1
erf.if_use_most = true
erf.if_z0 = 0.01
erf.use_coriolis = false

erf.ibseb.enable = true
erf.ibseb.prognostic = false
erf.ibseb.debug = false
erf.ibseb.csv_int = -1
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

# The ground: the force-restore balance on the sweep's shortwave, starting at erf.rad_t_sfc,
# with a thin skin (a small heat capacity), so the open ground warms by about a kelvin a
# step and a mean that counted the held footprints would show it.
erf.radiation.seb_prognostic_enable = true
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_surface_heat_capacity = 100.0
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1  = 1
erf.plot2d_vars_1 = seb_sw_net ibseb_ground_shadow seb_t_sfc
erf.radiation.diag_enable = true
erf.radiation.diag_stdout_enable = false
erf.radiation.diag_csv_enable = true
erf.radiation.diag_file = "radiation_diag.csv"
