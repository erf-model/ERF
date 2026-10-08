# The immersed-boundary surface energy balance on two levels: a 40 m cube
# (x 140-180 m, y 220-260 m) and a 60 m tower 80 m east of it (x 260-300 m) on
# a 480 m periodic domain, 20 m cells on level 0 and 10 m on level 1, under a
# low sun from the east (zenith 70 deg), so the tower shades the cube's
# east-facing wall. Level 1 holds the cube only (in_box from
# Exec/CanonicalTests/SEB/ibseb_refinement_box.py --region 120 200 200 280 --fit tight),
# so its rays must find the tower in the column map of level 0. The checker
# compares the cube's faces on level 1 against a run whose level 1 holds both
# buildings (erf.city.in_box_hi = 360 300).
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
amr.ref_ratio_vect = 2 2 1
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
erf.ibseb.lw_mode = "gray"
erf.ibseb.sky_emissivity = 0.83
