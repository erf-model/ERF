# Two-way coupled, fully compressible two-level ABL with many acoustic substeps (issue #5003).
#
# The divergence damping of the acoustic substeps has a diffusivity beta c^2 dtau, so with a
# fixed coefficient the damping per unit time falls as the number of substeps rises, and the
# acoustic waves inside the patch grow once the coarse level takes more than about 60 substeps.
# The coefficient is erf.beta_d divided by the horizontal acoustic Courant number of the
# substep.  At the 144 substeps fixed here, step 40 has max |w| = 0.004 m/s and a pressure
# perturbation of 0.02 Pa; with the unscaled coefficient they are 0.084 m/s and 43 Pa.
# The test requires |w| < 0.02 m/s at step 40 (add_test_field_bounds).
erf.prob_name = "ABL"

max_step = 40

amrex.fpe_trap_invalid = 1

geometry.prob_extent = 640 640 320
amr.n_cell           = 32  32  32
geometry.is_periodic = 1 1 0
amr.max_grid_size_x  = 16
amr.max_grid_size_y  = 16
amr.max_grid_size_z  = 32

amr.max_level = 1
amr.ref_ratio_vect = 2 2 1
amr.n_error_buf = 0
amr.blocking_factor = 2
amr.refine_whole_domain_dir = 2
erf.refinement_indicators = city
erf.city.max_level = 1
erf.city.in_box_lo = 160.0 160.0
erf.city.in_box_hi = 480.0 480.0

erf.fixed_dt = 1.0
erf.fixed_mri_dt_ratio = 144

zlo.type = "surface_layer"
zhi.type = "SlipWall"
erf.most.z0 = 0.1
erf.most.pot_surf_temp = 300.0

erf.init_type = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.use_gravity = true
erf.molec_diff_type = "None"
erf.les_type = "Smagorinsky"
erf.Cs = 0.17
erf.fix_random_seed = true

erf.use_coriolis = true
erf.latitude = 40.0
erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind = 5.0 0.0 0.0

erf.v = 0
erf.sum_interval = -1
erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1 = 40
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta pressure
