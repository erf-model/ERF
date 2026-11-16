# Motivation: a wrfinput file carries land and soil fields that only a land model used to
# ask for (LAI, VEGFRA, TMN, TSLB, SMOIS) and others that were never read (ZNT, ALBEDO,
# SEAICE, SNOW, SNOWH). With the two-stream force-restore balance and no land model the
# reader now takes them, and ALBBCK, each one optional. The file is written by
# make_wrfinput.py; every field varies in both horizontal directions and with the soil
# layer. check_log.py compares the component count, minimum and maximum the reader reports
# with the file.
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

# The force-restore balance without a land model: this is what asks for the fields
erf.radiation_model = "TwoStream"
erf.radiation.sw_enabled = true
erf.radiation.lw_enabled = true
erf.rad_t_sfc = 290.0
erf.radiation.seb_enable = true
erf.radiation.seb_prognostic_enable = true
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_t_deep_default = 290.0
