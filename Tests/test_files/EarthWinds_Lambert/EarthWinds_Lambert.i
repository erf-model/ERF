# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Earth-relative winds on a rotated (Lambert conformal) grid, from a wrfinput.
#
# The 25 km grid (TRUELAT1/2 = 30/60, STAND_LON = -120) turns by up to 11.4
# degrees from east at its corners. The run makes one step and writes:
#   - stations NW and NE placed by lat/long and SW placed by x/y, each with
#     x_velocity, y_velocity and the earth-relative u_east, v_north,
#     wind_speed and wind_dir (Output_Stations/*.dat);
#   - a 2D plotfile with lon_m and the sampled u_east and v_north at k = 10,
#     beside the 3D x_velocity and y_velocity.
# check_station_earth_winds.py and check_2d_earth_winds.py compare them with
# the Lambert rotation alpha = -n (lon - STAND_LON), n the cone factor.
# See README.md for the data this needs.
max_step = 1
start_datetime = "2024-08-05 00:00:00"

geometry.prob_extent = 2500000.0 2500000.0 10000.0
amr.n_cell           = 100 100 79
geometry.is_periodic = 0 0 0
amr.max_level        = 0
amr.max_grid_size    = 32

xlo.type = "Outflow"
xhi.type = "Outflow"
ylo.type = "Outflow"
yhi.type = "Outflow"
zlo.type = "surface_layer"
zhi.type = "SlipWall"
erf.most.z0   = 0.1
erf.most.zref = 10.0

erf.prob_name      = "WPS"
erf.init_type      = "WRFInput"
erf.nc_init_file_0 = "wrfinput_d01"
erf.use_real_bcs   = false
erf.terrain_type   = StaticFittedMesh
erf.moisture_model = MoistNoCondensation
erf.buoyancy_type  = 1
erf.les_type       = "Smagorinsky2D"
erf.Cs             = 0.1
erf.pbl_type       = "None"
erf.fixed_dt           = 5.0
erf.fixed_mri_dt_ratio = 6
erf.check_int          = -1

erf.plot_file_1 = plt
erf.plot_int_1  = 1
erf.plot_vars_1 = x_velocity y_velocity
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1  = 1
erf.plot2d_vars_1 = lon_m
erf.plot2d_level_sets_1 = k10
erf.plot2d.level_set.k10.coordinate = model_index
erf.plot2d.level_set.k10.values     = 10
erf.plot2d.level_set.k10.fields     = u_east v_north

erf.station_names = NW NE SW
erf.NW.field      = x_velocity y_velocity u_east v_north wind_speed wind_dir
erf.NW.lat        = 41.3960
erf.NW.long       = -132.2761
erf.NW.height_agl = 100.0 1500.0
erf.NE.field      = x_velocity y_velocity u_east v_north wind_speed wind_dir
erf.NE.lat        = 39.4627
erf.NE.long       = -109.3217
erf.NE.height_agl = 100.0 1500.0
erf.SW.field      = x_velocity y_velocity u_east v_north wind_speed wind_dir
erf.SW.x          = 262500.0
erf.SW.y          = 262500.0
erf.SW.height_agl = 100.0 1500.0
erf.station_sampling_interval = 1
