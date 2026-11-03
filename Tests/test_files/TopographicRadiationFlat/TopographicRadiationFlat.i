# Flat fitted terrain makes the slope factor exactly one and casts no shadows.
# Compare the active correction with the ordinary RRTMGP-to-SLM path after
# both the first radiation call and the next scheduled update.
erf.prob_name = ABL
max_step = 3
stop_time = 1.0e6
# HDF5 initializes its native float types with invalid-operation probes.
amrex.fpe_trap_invalid = 0

geometry.prob_lo = 0 0 0
geometry.prob_hi = 400 400 16000
amr.n_cell = 8 8 32
geometry.is_periodic = 1 1 0
amr.max_level = 0

erf.terrain_type = StaticFittedMesh
erf.flat_terrain = true
zlo.type = surface_layer
zhi.type = SlipWall
erf.surface_layer.flux_type = custom
erf.most.z0 = 0.1
erf.most.zref = 25.0
erf.most.ustar = 0.1
erf.most.tstar = 0.0
erf.most.qstar = 0.0

erf.moisture_model = Kessler
erf.land_surface_model = SLM
slm.nsoil = 4
slm.soil_dz = 0.05 0.10 0.30 0.55
slm.landtype0 = 10
slm.LAI0 = 2.0
slm.clay0 = 13.0
slm.sand0 = 17.0
slm.sw0 = 0.60 0.65 0.70 0.80
slm.st0 = 300.0 299.8 299.5 299.0
slm.relax_hgt = 0.0 0.0 0.0 1.0
slm.soiltnudging = true
slm.soilwnudging = true
slm.tausoil = 86400.0
slm.tabs_s = 0.0
slm.t00 = 300.0
slm.z0_soil = 0.0387
slm.mws_mx0 = 50.0
slm.Rc_max = 5000.0
slm.T_opt = 298.0

erf.init_type = input_sounding
erf.input_sounding_file = input_sounding
erf.use_gravity = true
erf.fixed_dt = 0.5
erf.les_type = Smagorinsky
erf.Cs = 0.1
erf.molec_diff_type = None

start_datetime = "2024-08-05 18:00:00"
erf.radiation_model = RRTMGP
erf.rad_t_sfc = 300.0
erf.rad_cons_lat = 40.0
erf.rad_cons_lon = -100.0
erf.rad_freq_in_steps = 2
erf.rad_do_subcol_sampling = false
erf.rrtmgp_file_path = radiation_data
erf.rrtmgp_coeffs_sw = rrtmgp-data-sw-g224-2018-12-04.nc
erf.rrtmgp_coeffs_lw = rrtmgp-data-lw-g256-2018-12-04.nc
erf.rrtmgp_cloud_optics_sw = rrtmgp-cloud-optics-coeffs-sw.nc
erf.rrtmgp_cloud_optics_lw = rrtmgp-cloud-optics-coeffs-lw.nc
erf.co2vmr = 0.00042
erf.o3vmr = 3.0e-8
erf.n2ovmr = 3.2e-7
erf.covmr = 1.5e-7
erf.ch4vmr = 1.7e-6
erf.o2vmr = 0.209
erf.n2vmr = 0.7906
erf.fixed_solar_zenith_angle = 0.5
erf.fixed_solar_azimuth_angle = 90

erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1 = 1
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta qv qsrc_sw qsrc_lw
erf.plot_lsm = true
erf.plot2d_file_2 = terrain
erf.plot2d_int_2 = 1
erf.plot2d_vars_2 = topographic_shadow terrain_direct_factor unshaded_sw_flux_dn terrain_sw_flux_dn
