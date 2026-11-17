# WRFInput_LegacySurfaceTemperature.i with a level 1 that lists two init files (see
# make_wrfinputs.py): nest_d02, read at the start, and nest_d03, which starts 10 s later
# and so is never read by a run as short as this.  The initial run writes nest_chk00000.
max_step = 0
start_datetime = "2026-01-01 00:00:00"
geometry.prob_extent = 16000.0 16000.0 800.0
amr.n_cell = 16 16 8
geometry.is_periodic = 0 0 0
amr.max_level = 1
amr.ref_ratio_vect = 3 3 1
amr.max_grid_size = 64
amr.refine_whole_domain_dir = 2

xlo.type = "Outflow"
xhi.type = "Outflow"
ylo.type = "Outflow"
yhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.prob_name = "WPS"
erf.init_type = "WRFInput"
erf.nc_init_file_0 = "nest_d01"
erf.nc_init_file_1 = "nest_d02" "nest_d03"
erf.use_real_bcs = false
erf.terrain_type = StaticFittedMesh
erf.avg_grid_faces_to_nodes = true
erf.rebalance_wrf_input = false
erf.moisture_model = MoistNoCondensation
erf.pbl_type = None
erf.les_type = None
erf.coupling_type = OneWay
erf.fixed_dt = 1.0
erf.check_file = nest_chk
erf.check_int = 1
