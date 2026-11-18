# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# The two-stream force-restore balance is a model of soil. It used to run on every
# column: on water (land mask 0) the skin and the soil water evolved as soil, and the
# sweep took that skin as its longwave boundary, while the surface layer kept the sea
# surface temperature there. The balance now reads the land mask: a water column keeps
# its state, is left out of the means, and the sweep takes the surface layer's sea
# surface temperature.
#
# RunTwoStreamSEBLandMask.cmake runs this deck with erf.is_land = 1 and = 0
# (an idealized run has one mask value per level), and a third time with two levels,
# water on level 0 and a level 1 of land over its western half. check_land_mask.py asserts
#   1. water: seb_t_sfc and seb_q_sfc stay at their starting values at every step;
#   2. land: seb_t_sfc moves (so check 1 is not a run in which nothing evolves);
#   3. step 1: the net longwave of the two runs differs by
#      sigma (T_sea^4 - T_skin^4), with T_sea = t_surf (surf_pres / p0)^(R/cp) from the
#      water run's own 2D output and T_skin = erf.rad_t_sfc;
#   4. the radiation diagnostics file reports a mean skin on land and none on water;
#   5. two levels: the level 0 water columns under the level 1 land keep their state
#      (the average-down between the levels is over land columns), and level 1 evolves.
#
# This is the deck of TwoStream_SEBSurfaceLayerFluxes with a surface emissivity of one,
# so that check 3 has no reflected term. The surface pressure is 950 hPa and
# erf.most.pot_surf_temp is 301.5 K, so the sea is at about 297.1 K against a skin
# of 299 K: the two longwave boundaries differ by about 11 W/m^2.

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
erf.plot2d_vars_1 = seb_t_sfc seb_q_sfc seb_lw_net t_surf surf_pres

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
# fluxes. The feature under test is where it runs: on land columns only.
erf.radiation.seb_enable = true
erf.radiation.seb_prognostic_enable = true
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_turbulent_flux_source = surface_layer
erf.radiation.seb_t_deep_default = 300.0
erf.radiation.seb_q_sfc_default = 0.0095
erf.radiation.seb_q_deep_default = 0.0095
erf.radiation.seb_surface_heat_capacity = 2.0e4

erf.radiation.surface_emissivity_lw = 1.0
erf.radiation.diag_enable = true
erf.radiation.diag_file = "radiation_diag.dat"
