# REQ-SHOC-01/02/03: explicit native SHOC state-update ownership.
# The finite near-saturation layer in input_sounding_stable_cloud is a
# deliberately preconditioned idealized cloud fixture, not a BOMEX contract.
erf.prob_name = "ABL"
erf.init_type = "input_sounding"
erf.sounding_type = Ideal
erf.input_sounding_file = "input_sounding_stable_cloud"
erf.moisture_model = SatAdj

stop_time = 20.0
max_step = 20
erf.fixed_dt = 1.0
erf.fixed_mri_dt_ratio = 2
erf.fix_random_seed = 1

geometry.prob_extent = 400.0 400.0 480.0
amr.n_cell = 4 4 48
amr.max_level = 0
amr.max_grid_size = 64
fabarray.mfiter_tile_size = 1024 1024 1024
geometry.is_periodic = 1 1 0
zlo.type = "surface_layer"
zhi.type = "SlipWall"

erf.surface_layer.flux_type = "custom"
erf.most.ustar = 0.25
erf.most.tstar = 0.005
erf.most.qstar = -0.00002
erf.most.z0 = 0.10
erf.most.zref = 10.0

erf.molec_diff_type = "None"
erf.rans_type = "None"
erf.les_type = "Smagorinsky2D"
erf.Cs = 0.4
erf.mix_isotropic = false
erf.use_Ri_correction = true
erf.pbl_type = "NATIVE_SHOC"

erf.shoc.transport_mode = state_update
erf.shoc.momentum_transport = state_update
erf.shoc.signed_tke_production = true
erf.shoc.extra_shoc_diags = true
erf.shoc.shoc_1p5tke = true
erf.shoc.allow_tendency_microphysics_overlap = false
erf.shoc.top_taper_depth = 120.0
erf.shoc.top_taper_min_factor = 0.25
erf.shoc.debug_summary = false
erf.shoc.check_flux_state = true
erf.shoc.column_conservation_check = true
erf.shoc.debug_bad_column = false
erf.shoc.debug_bad_column_abort = false
erf.shoc.debug_disable_pdf_cloud_increment = false
erf.shoc.debug_disable_theta_state_update = false
erf.shoc.debug_disable_moisture_state_update = false
erf.shoc.debug_disable_tke_state_update = false

erf.vert_implicit = false
erf.use_gravity = true
erf.implicit_before_substep = false
erf.implicit_thermal_diffusion = false
erf.implicit_moisture_diffusion = false
erf.implicit_ke_diffusion = false
erf.implicit_momentum_diffusion = false
erf.vert_implicit_fac = 0.0 0.0 0.0

erf.plot_file_1 = plt
erf.plot_int_1 = 10
erf.plot_vars_1 = density temp theta pressure qv qc rhoKE x_velocity y_velocity Kmv Khv Lturb shoc_cldfrac shoc_ql shoc_cond brunt shear_prod buoy_prod diss_tke

amrex.fpe_trap_invalid = 1
erf.v = 0
amr.v = 0
