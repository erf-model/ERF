# Build an 8x8 legacy-checkpoint fixture. The restart migration tests attach a larger
# WRFInput file and verify that only the local PSFC region is recovered.
max_step = 0
start_datetime = "2026-01-01 00:00:00"
geometry.prob_extent = 8000.0 8000.0 800.0
amr.n_cell = 8 8 8
geometry.is_periodic = 0 0 0
amr.max_level = 0
amr.blocking_factor = 2
amr.max_grid_size = 4 4 8

xlo.type = "Outflow"
xhi.type = "Outflow"
ylo.type = "Outflow"
yhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.prob_name = "WPS"
erf.init_type = "WRFInput"
erf.nc_init_file_0 = "wrfinput_crop_small"
erf.use_real_bcs = false
erf.terrain_type = StaticFittedMesh
erf.avg_grid_faces_to_nodes = true
erf.rebalance_wrf_input = false
erf.moisture_model = MoistNoCondensation
erf.pbl_type = None
erf.les_type = None
erf.fixed_dt = 1.0
erf.check_file = crop_chk
erf.check_int = 1
