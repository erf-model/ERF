# A one-level WRFInput run whose wrfinput PSFC, SST and TSK vary in space (see
# make_wrfinputs.py).  The initial run writes chk00000; legacy_restarts.py turns copies
# of it into the checkpoints older versions wrote, and restarts from each.
max_step = 0
start_datetime = "2026-01-01 00:00:00"
geometry.prob_extent = 16000.0 16000.0 800.0
amr.n_cell = 16 16 8
geometry.is_periodic = 0 0 0
amr.max_level = 0
amr.max_grid_size = 8 8 8

xlo.type = "Outflow"
xhi.type = "Outflow"
ylo.type = "Outflow"
yhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.prob_name = "WPS"
erf.init_type = "WRFInput"
erf.nc_init_file_0 = "wrfinput_d01"
erf.use_real_bcs = false
erf.terrain_type = StaticFittedMesh
erf.avg_grid_faces_to_nodes = true
erf.rebalance_wrf_input = false
erf.moisture_model = MoistNoCondensation
erf.pbl_type = None
erf.les_type = None
erf.fixed_dt = 1.0
erf.check_file = chk
erf.check_int = 1
