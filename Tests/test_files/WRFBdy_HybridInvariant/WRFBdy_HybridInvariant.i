# Motivation: the wrfbdy arrays are the WRF state coupled with the column mass and the map
# factors -- PH with c1f mu + c2f on full levels (issue #5149), U and V with
# (c1h mu + c2h)/msfuy and /msfvx. At the first boundary time, once ERF has decoupled them
# and interpolated them onto its heights, they must give back the initial state that
# init_from_wrfinput builds from wrfinput_d01. make_wrf_hybrid.py writes such a pair on the
# hybrid coordinate (hybrid_opt = 2) with terrain and map factors other than one;
# erf.check_wrfbdy_invariant prints the largest difference on each side and
# check_invariant.py compares it with the tolerances.
max_step = 0
geometry.prob_extent = 48000.0 36000.0 20000.0
amr.n_cell = 16 12 20
geometry.is_periodic = 0 0 0
amr.max_level = 0
amr.max_grid_size = 8

zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.prob_name = "WPS"
erf.init_type = "WRFInput"
erf.nc_init_file_0 = "wrfinput_d01"
erf.nc_bdy_file = "wrfbdy_d01"
erf.real_width = 3
erf.terrain_type = StaticFittedMesh
# The ERF heights come from the WRF ones and the state is not rebalanced, so that the
# initial density is 1/ALT, as the boundary density from the geopotential is
erf.avg_grid_faces_to_nodes = true
erf.rebalance_wrf_input = false
erf.moisture_model = "Kessler"
erf.pbl_type = None
erf.les_type = None
erf.fixed_dt = 1.0
erf.check_int = -1
erf.plot_int_1 = -1

# What this case tests
erf.check_wrfbdy_invariant = true
