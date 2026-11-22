# Motivation: the two-stream force-restore balance started its skin temperature at the one
# value erf.rad_t_sfc in every column. With erf.radiation.seb_init_skin_from_file a land
# column starts at the TSK of the wrfinput file and a water column keeps erf.rad_t_sfc.
# make_wrfinput.py writes TSK = 280 + 0.5 i + 0.25 j [K] and a block of water columns that
# lies in both boxes; erf.rad_t_sfc is 300 K, a value TSK never takes. check_skin.py
# compares the skin of every column of the 2D plotfile with the formula.
max_step = 0
start_datetime = "2026-01-01 00:00:00"
geometry.prob_extent = 12000.0 8000.0 800.0
amr.n_cell = 12 8 8
geometry.is_periodic = 0 0 0
amr.max_level = 0
amr.max_grid_size = 8

xlo.type = "Outflow"
xhi.type = "Outflow"
ylo.type = "Outflow"
yhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.prob_name = "WPS"
erf.init_type = "WRFInput"
erf.nc_init_file_0 = "wrfinput_d01"
erf.nc_low_file = "wrflowinp_d01"
erf.use_real_bcs = false
erf.terrain_type = StaticFittedMesh
erf.avg_grid_faces_to_nodes = true
erf.rebalance_wrf_input = false
erf.moisture_model = MoistNoCondensation
erf.pbl_type = None
erf.les_type = None
erf.fixed_dt = 1.0
erf.check_int = -1
erf.plot_int_1 = -1

# The skin and the deep temperature of every column, at the start
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 1
erf.plot2d_vars_1 = seb_t_sfc seb_t_deep

# The force-restore balance without a land model
erf.radiation_model = "TwoStream"
erf.radiation.sw_enabled = true
erf.radiation.lw_enabled = true
erf.rad_t_sfc = 300.0
erf.radiation.seb_enable = true
erf.radiation.seb_prognostic_enable = true
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_t_deep_default = 290.0

# What this case tests
erf.radiation.seb_init_skin_from_file = true
# The WRFInput_SEBDeepInit legs set this one: a land column's deep temperature is then the
# TMN of the file, 275 + 0.25 i + 0.5 j [K], and a water column keeps seb_t_deep_default
erf.radiation.seb_init_t_deep_from_file = false
