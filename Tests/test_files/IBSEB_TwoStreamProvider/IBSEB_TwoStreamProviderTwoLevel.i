# Two levels on the two-stream columns, refined in x, y AND z: the deck of
# Tests/test_files/IBSEB_RefinedLevels (a 40 m cube on level 1, a 60 m tower outside
# it, 20 m cells on level 0) with level 1 refined by 2 in every direction (10 m by
# 10 m by 5 m cells, 32 layers against 16), the two-stream radiation, the faces'
# radiation from its columns (erf.ibseb.radiation = two_stream) and their sky and
# ground longwave from them too (erf.ibseb.lw_mode = two_stream, the default with that
# provider). The sun is fixed 70 degrees from the zenith, in the east, for the faces
# and the columns alike. Tests/RunIBSEBTwoStreamProvider.cmake runs it under a
# transparent and an absorbing sky.
max_step  = 0
amrex.fpe_trap_invalid = 1
amrex.fpe_trap_zero = 1
amrex.fpe_trap_overflow = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_extent =  480     480    160
amr.n_cell           =  24      24     16
amr.max_grid_size    =  64
geometry.is_periodic = 1 1 0
zlo.type = "surface_layer"
erf.most.z0 = 0.1
erf.most.surf_temp = 300.0
zhi.type = "SlipWall"
erf.fixed_dt           = 0.5
erf.substepping_cfl    = 0.5
erf.sum_interval   = -1
erf.v              = 0
amr.v              = 0

amr.max_level = 1
amr.ref_ratio_vect = 2 2 2
amr.n_error_buf = 0
amr.blocking_factor = 2
amr.refine_whole_domain_dir = 2
erf.refinement_indicators = city
erf.city.max_level = 1
erf.city.in_box_lo = 100.0 180.0
erf.city.in_box_hi = 220.0 300.0

erf.check_file      = chk
erf.check_int       = -1
erf.plot_file_1     = plt
erf.plot_int_1      = -1
erf.use_gravity = true
erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.17
erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.buildings_type = ImmersedForcing
erf.buildings_file_name = cube_and_tower_10m.txt
erf.immersed_forcing_substep = true
eb2.small_volfrac = 0.005
erf.if_use_most = true
erf.if_z0 = 0.01
erf.if_snap_partial_cells = true
erf.use_coriolis = false

erf.ibseb.enable = true
erf.ibseb.prognostic = true
erf.ibseb.csv_int = 1
erf.ibseb.csv_file = "ibseb_buildings.csv"
erf.ibseb.dump_faces_file = "faces/set"
erf.ibseb.T_skin_init = 300.0
erf.ibseb.T_interior = 295.0
erf.ibseb.n_slab_layers = 4
erf.ibseb.view_n_az = 16
erf.ibseb.view_n_el = 8
erf.ibseb.sun_mode = "fixed"
erf.ibseb.sun_zenith_deg = 70.0
erf.ibseb.sun_azimuth_deg = 90.0
erf.ibseb.radiation = "two_stream"

erf.radiation_model = TwoStream
erf.fixed_solar_zenith_angle = 0.3420201433256688
erf.fixed_total_solar_irradiance = 1000.0
erf.rad_t_sfc = 300.0
