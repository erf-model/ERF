# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# erf.radiation.seb_snow_enable puts one bulk snow layer on the two-stream
# force-restore surface (Source/Radiation/TwoStream/ERF_ForceRestoreSnow.H), the
# snow routine of the ISBA scheme as in Douville, Royer and Mahfouf (1995,
# Climate Dynamics 12, 21-35), Sec. 2.1, for bare ground.
#
# Tests/RunTwoStreamSEBSnow.cmake runs this deck in several legs and
# check_two_stream_seb_snow.py (this folder) checks them against the paper's
# closed forms; README.md lists the legs and the checks.
#
# The column is the TwoStream_SEBSurfaceLayerFluxes one: Kessler, a fixed sun,
# surface-layer fluxes, and a skin that starts at 299 K, above the melting
# point, so a snowpack given at the start melts.

erf.prob_name = "ABL"

max_step = 10
stop_time = 1.0e6
amrex.fpe_trap_invalid = 0

geometry.prob_extent = 400 400 800
amr.n_cell           = 4 4 16
geometry.is_periodic = 1 1 0

zlo.type = "surface_layer"
erf.most.z0        = 0.1
erf.most.zref      = 25.0
erf.most.pot_surf_temp = 301.5     # warmer than the air: H > 0
erf.most.surf_moist = 0.0095   # moister than the air: LE > 0
zhi.type = "SlipWall"
zhi.theta_grad = 0.003

erf.fixed_dt = 1.0
erf.v = 0
amr.v = 0
amr.max_level = 0

erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1 = -1
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 1
erf.plot2d_vars_1 = seb_t_sfc seb_q_sfc seb_lh seb_sw_net seb_snow_swe seb_snowfall_total seb_snow_sublimation_total seb_snow_melt_total seb_rain_total seb_runoff_total precip_rain_accum precip_frozen_accum

erf.use_gravity = true
erf.molec_diff_type = "None"
erf.les_type = "None"
erf.pbl_type = "MRF"
erf.theta_ref = 300.0
erf.moisture_model = "Kessler"

erf.init_type = "input_sounding"
erf.sounding_type = Ideal
erf.input_sounding_file = "input_sounding"
erf.use_coriolis = false
erf.abl_driver_type = "None"

# RADIATION - TwoStream, SW + LW, clear sky, fixed sun
erf.radiation_model = "TwoStream"
erf.radiation.sw_enabled = true
erf.radiation.lw_enabled = true
erf.radiation.tau_per_layer = 0.00625
erf.radiation.tau_lw_per_layer = 1.0
erf.fixed_solar_zenith_angle = 0.5    # cos(60 deg)
erf.fixed_total_solar_irradiance = 1361.0
erf.rad_t_sfc = 299.0    # the skin's starting temperature [K]
erf.radiation.v = 0

# The prognostic surface energy balance, driven by the sweep's own surface
# fluxes and the surface layer's H and LE.
erf.radiation.seb_enable = true
erf.radiation.seb_prognostic_enable = true
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_turbulent_flux_source = surface_layer
erf.radiation.seb_t_deep_default = 300.0
erf.radiation.seb_q_sfc_default = 0.0095
erf.radiation.seb_q_deep_default = 0.0095
erf.radiation.seb_surface_heat_capacity = 2.0e4

# The force-restore hydrology takes the rain and the snow's melt water.
erf.radiation.seb_precipitation_source = microphysics
erf.radiation.seb_soil_porosity = 0.45

# The snow layer under test (Douville et al. 1995: W_crn = 10 kg/m^2, albedo 0.70).
erf.radiation.seb_snow_enable = true
erf.radiation.seb_snow_initial_swe_kg_m2 = 50.0
erf.radiation.seb_snow_albedo = 0.70
erf.radiation.seb_snow_critical_swe_kg_m2 = 10.0
