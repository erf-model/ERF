# Lines given by latitude and longitude (erf.sample_line_lo/hi_latlon), on the rotated
# lat/lon grid of the wrfinput that make_wrfinput.py writes. theta = 300 + 0.1 i + 0.001 j
# names the column each sample comes from, and check_line_latlon.py requires:
#   - prof: a vertical line at the centre of cell (10,20), all 8 cells of the column;
#   - ew:   along x from the centre of (4,5) to a point 0.4 cells north of the centre of
#           (25,5), sampled along row j = 5 (i = 4..25);
#   - ns:   along y from the centre of (7,3) to a point 0.7 cells west of the centre of
#           (7,28), sampled along column i = 7 (j = 3..28).
# erf.sample_line_latlon places the column of the log col_ll at the same site as prof, which
# must be the centre of cell (10,20). The boxes are 16 x 16, so each horizontal line crosses
# two of them.
max_step = 1
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

erf.do_line_sampling          = true
erf.line_sampling_interval    = 1
erf.line_sampling_text_output = true
erf.line_sampling_vars        = theta
erf.sample_line_name          = prof ew ns
erf.sample_line_dir           = 2 0 1
erf.sample_line_lo_latlon     = 40.205500 -104.918150  50.   40.058500 -104.958350 250.   40.046500 -104.919250 250.
erf.sample_line_hi_latlon     = 40.205500 -104.918150 750.   40.104100 -104.713450 250.   40.270100 -104.977440 250.

erf.sum_interval       = 1
erf.sample_line_log    = col_ll
erf.sample_line_latlon = 40.205500 -104.918150
