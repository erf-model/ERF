# erf.sample_point and its forms (#4978), on the wrfinput that make_wrfinput.py --hill writes:
# a rotated lat/lon grid whose ground rises 20 m per cell in i, under a flat top at 800 m, with
# 8 levels spaced evenly in each column, theta naming each column and qv = 0.001 (k+1) each
# level. pt_a is in column (10,20), whose ground is at 210 m and cells 73.75 m deep; pt_b in
# column (25,5), ground at 510 m and cells 36.25 m deep. The test gives the points as
#   erf.sample_point        = 10 20 3   25 5 6                               (k = 3, 6)
#   erf.sample_point_real   = 10500 20500 560   25500 5500 700               (k = 4, 5)
#   erf.sample_point_real   = 10500 20500 50    25500 5500 200, above ground (k = 0, 5)
#   erf.sample_point_latlon = the two cell centres, 250 and 200 above ground (k = 3, 5)
# check_sample_point.py requires each point's log to hold the conserved state of the expected
# cell of the full-depth columns col_a and col_b, placed by erf.sample_line_latlon at the same
# cell centres. The boxes are 16 x 16, so the points are on different boxes.
max_step = 2
start_datetime = "2026-01-01 00:00:00"
geometry.prob_extent = 32000.0 32000.0 800.0
amr.n_cell           = 32 32 8
geometry.is_periodic = 0 0 0
amr.max_level        = 0
amr.max_grid_size    = 16

xlo.type = "Outflow"
xhi.type = "Outflow"
ylo.type = "Outflow"
yhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.prob_name       = "WPS"
erf.init_type       = "WRFInput"
erf.nc_init_file_0  = "wrfinput_d01"
erf.nc_low_file     = "wrflowinp_d01"
erf.use_real_bcs    = false
erf.terrain_type    = StaticFittedMesh
erf.moisture_model  = MoistNoCondensation
erf.pbl_type        = None
erf.les_type        = None
erf.fixed_dt        = 1.0
erf.check_int       = -1

erf.sum_interval       = 1
erf.sample_line_log    = col_a col_b
erf.sample_line_latlon = 40.205500 -104.918150   40.100500 -104.712650
erf.sample_point_log   = pt_a pt_b
