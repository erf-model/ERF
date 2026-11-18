# d02's full footprint begins at level-1 index 18; the ERF patch begins at 24.
# d03 must begin at level-2 index 78, independent of the level-1 crop.
max_step = 0
start_datetime = "2026-01-01 00:00:00"
geometry.prob_extent = 48000.0 48000.0 800.0
amr.n_cell = 48 48 8
geometry.is_periodic = 0 0 0
amr.max_level = 2
amr.ref_ratio_vect = 3 3 1 3 3 1
amr.max_grid_size = 64
amr.refine_whole_domain_dir = 2

erf.refinement_indicators = crop
erf.crop.max_level = 1
erf.crop.in_box_lo_indices = 24 24 0
erf.crop.in_box_hi_indices = 47 47 7

xlo.type = "Outflow"
xhi.type = "Outflow"
ylo.type = "Outflow"
yhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.prob_name = "WPS"
erf.init_type = "WRFInput"
erf.nc_init_file_0 = "wrfinput_d01"
erf.nc_init_file_1 = "wrfinput_d02"
erf.nc_init_file_2 = "wrfinput_d03"
erf.nc_low_file = "wrflowinp_d01"
erf.nc_low_file_1 = "wrflowinp_d02"
erf.nc_low_file_2 = "wrflowinp_d03"
erf.use_real_bcs = false
erf.terrain_type = StaticFittedMesh
erf.avg_grid_faces_to_nodes = true
erf.rebalance_wrf_input = false
erf.moisture_model = MoistNoCondensation
erf.pbl_type = None
erf.les_type = None
erf.coupling_type = OneWay
erf.fixed_dt = 1.0
erf.check_file = chk
erf.check_int = 1

erf.plot_file_1 = plt
erf.plot_int_1 = 1
erf.plot_vars_1 = theta
