# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# TwoStream_SEBSurfaceLayerFluxes on a refined hierarchy. Each level's surface
# energy balance takes the H and LE its own surface layer applies (the level's
# SFS fluxes), and with erf.radiation.seb_surface_layer_uses_skin each level's
# surface layer takes that level's skin. The single-level deck cannot see a level
# reading another level's fields, or a level created mid-run starting without them.
#
# A Witch-of-Agnesi ridge makes every column different, so a flux or skin taken
# from the wrong column or level shows. Level 1 covers the middle half of the
# domain in x, spans z (amr.refine_whole_domain_dir = 2, which the MRF column
# scheme needs), and is created by the regrid at the start of step 7, the first
# step that starts after erf.patch.start_time (erf.fixed_dt = 0.5 s).
#
# check_two_stream_seb_flux_source.py --multilevel asserts, on every level and
# every step, that seb_hfx / seb_lh equal sensible_heat_flux / latent_heat_flux
# column by column, and in the two-way run that t_surf is the previous step's
# skin as a potential temperature. Its budget check is single-level only: under
# the patch the coarse skin is the average of the fine one.

erf.prob_name = "ABL"

max_step = 10
stop_time = 1.0e6
amrex.fpe_trap_invalid = 0

geometry.prob_extent = 800 400 800
amr.n_cell           = 8 4 16
geometry.is_periodic = 1 1 0

zlo.type = "surface_layer"
erf.most.z0        = 0.1
erf.most.zref      = 25.0
erf.most.average_policy = 1   # per-column u*, theta*: H differs from column to column
erf.most.surf_temp = 301.5     # warmer than the air: H > 0
erf.most.surf_moist = 0.0095   # moister than the air: LE > 0
zhi.type = "SlipWall"
zhi.theta_grad = 0.003

erf.fixed_dt = 0.5     # the CTest runner passes the same value (DT 0.5)
erf.v = 0
amr.v = 0
amr.max_level = 1
amr.ref_ratio_vect = 2 2 1
amr.n_error_buf = 0
amr.blocking_factor = 2
amr.refine_whole_domain_dir = 2
erf.refinement_indicators = patch
erf.patch.max_level = 1
erf.patch.in_box_lo = 200.0 0.0
erf.patch.in_box_hi = 600.0 400.0
erf.patch.start_time = 2.75    # level 1 appears at step 7, which starts at t = 3.0 (dt 0.5 s)
erf.regrid_int = 1
# No subcycling: with the default two fine substeps per coarse step, the fine surface
# layer's t_surf at an output comes from the skin of the last fine substep, not from the
# skin the previous output holds, and the two-way check could not pair them.
erf.dt_ref_ratio = 1

# A ridge across x: every column has its own height and lowest-cell thickness.
erf.terrain_type         = StaticFittedMesh
erf.terrain_smoothing    = 0
prob.custom_terrain_type = "WoA"
prob.dir                 = 0
prob.hmax                = 40.0
prob.L                   = 150.0

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
