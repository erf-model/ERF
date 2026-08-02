erf.prob_name = "CustomSurfaceFluxFile"

max_step = 2

geometry.prob_extent = 80 80 400
amr.n_cell = 8 8 16
geometry.is_periodic = 1 1 0

zlo.type = "surface_layer"
zhi.type = "SlipWall"
erf.surface_layer.flux_type = "custom_file"
erf.most.flux_file_0 = "surface_flux_0.nc"
erf.most.flux_file_1 = "surface_flux_1.nc"
erf.most.zref = 10.0

erf.init_type = "input_sounding"
erf.input_sounding_file = "input_sounding"

erf.fixed_mri_dt_ratio = 2
erf.fixed_dt = 1.0
erf.dt_ref_ratio = 1
erf.sum_interval = 1
erf.v = 1
amr.v = 1

amr.max_level = 1
amr.ref_ratio_vect = 2 2 1
erf.coupling_type = "OneWay"
erf.refinement_indicators = box1
erf.box1.max_level = 1
erf.box1.in_box_lo = 20.0 20.0 0.0
erf.box1.in_box_hi = 60.0 60.0 400.0
amr.n_error_buf = 0
amr.blocking_factor = 1

erf.check_int = -1
erf.plot_int_1 = -1
erf.plotfile2d_type = "Netcdf"
erf.plot2d_file_1 = "surface_"
erf.plot2d_int_1 = 2
erf.plot2d_vars_1 = u_star t_star q_star sens_flux laten_flux

erf.moisture_model = "SatAdj"
erf.vert_implicit_fac = 0.0 0.0 0.0
erf.molec_diff_type = "Constant"
erf.alpha_T = 1.0
erf.alpha_C = 1.0
erf.dynamic_viscosity = 1.0
erf.use_gravity = true
erf.use_coriolis = false
erf.les_type = "None"
erf.pbl_type = "None"
