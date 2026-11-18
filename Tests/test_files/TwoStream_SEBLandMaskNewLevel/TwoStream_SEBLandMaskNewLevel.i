max_step = 8
start_datetime = "2026-06-21 12:00:00"
geometry.prob_extent = 16000.0 16000.0 800.0
amr.n_cell           = 16 16 8
geometry.is_periodic = 0 0 0
amr.max_grid_size    = 8

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
erf.v = 0
amr.v = 1

# Level 1 over coarse columns i = 6..11, j = 4..9, tagged from t = 4.5 s on: the regrid
# after step 5 makes it (ERF::MakeNewLevelFromCoarse). The box is given in level 1 indices.
amr.max_level = 1
amr.ref_ratio_vect = 2 2 1
amr.n_error_buf = 0
amr.blocking_factor = 4
amr.refine_whole_domain_dir = 2
erf.regrid_int = 1
erf.dt_ref_ratio = 1
# Two-way: the average-down of the surface fields runs between the levels at every step.
erf.coupling_type = TwoWay
erf.refinement_indicators = patch
erf.patch.max_level = 1
erf.patch.in_box_lo_indices = 12 8 0
erf.patch.in_box_hi_indices = 23 19 7
erf.patch.start_time = 4.5

erf.check_file = chk
erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1 = -1
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 1
erf.plot2d_vars_1 = seb_t_sfc seb_q_sfc

erf.radiation_model = "TwoStream"
erf.radiation.sw_enabled = true
erf.radiation.lw_enabled = true
erf.radiation.tau_per_layer = 0.00625
erf.radiation.tau_lw_per_layer = 1.0
erf.fixed_solar_zenith_angle = 0.5
erf.fixed_total_solar_irradiance = 1361.0
erf.rad_t_sfc = 299.0
erf.radiation.v = 0

erf.radiation.seb_enable = true
erf.radiation.seb_prognostic_enable = true
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_turbulent_flux_source = defaults
erf.radiation.seb_t_deep_default = 300.0
erf.radiation.seb_q_sfc_default = 0.0095
erf.radiation.seb_q_deep_default = 0.0095
erf.radiation.seb_surface_heat_capacity = 2.0e4
