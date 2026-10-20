# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# The prognostic surface energy balance removes the sensible (H) and latent (LE)
# heat fluxes from the ground. Without a land model it used to take them from
# erf.radiation.seb_hfx_default and seb_lh_default (0 by default), while the
# surface layer was putting its own H and LE into the air: the ground never lost
# what the air gained, and the skin ran hot by day. With
# erf.radiation.seb_turbulent_flux_source = surface_layer (the default) the
# balance takes the surface layer's applied fluxes.
#
# check_two_stream_seb_flux_source.py runs this deck three times: surface_layer
# (one-way), defaults, and surface_layer with
# erf.radiation.seb_surface_layer_uses_skin = true (two-way). It asserts
#   1. one-way and two-way: the balance's H and LE (seb_hfx, seb_lh) equal the
#      surface layer's sensible_heat_flux and latent_heat_flux to round-off at
#      every step, and are not small;
#   2. defaults: they are the constants (0);
#   3. one-way and two-way: the skin ends cooler than with the defaults by the
#      energy those fluxes carried away, sum(dt (H + LE)) / C_s, to within 5 %;
#   4. one-way: the surface layer keeps its own surface temperature
#      (erf.most.surf_temp, a potential temperature);
#   5. two-way: the surface layer's t_surf is the skin temperature of the step
#      before, converted to potential temperature, (p0 / p_sfc)^(R/cp) T_s.
#
# The surface pressure is 950 hPa so that conversion is a 1.5 % factor, and the
# skin starts at 299 K, away from the surface layer's 301.5 K, so a skin that is
# never handed over or not converted fails assertion 5.

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
erf.most.surf_temp = 301.5     # warmer than the air: H > 0
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
erf.plot2d_vars_1 = seb_t_sfc seb_hfx seb_lh sensible_heat_flux latent_heat_flux t_surf surf_pres

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
# fluxes. The feature under test is its source of H and LE.
erf.radiation.seb_enable = true
erf.radiation.seb_prognostic_enable = true
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_turbulent_flux_source = surface_layer
erf.radiation.seb_t_deep_default = 300.0
erf.radiation.seb_q_sfc_default = 0.0095
erf.radiation.seb_q_deep_default = 0.0095
erf.radiation.seb_surface_heat_capacity = 2.0e4
