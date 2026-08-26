# Minimal deterministic NetCDF terrain + gridded forest regression.
geometry.prob_lo     = 0. 0. 0.
geometry.prob_hi     = 2000. 2000. 1000.
amr.n_cell           = 16 16 16
amr.max_level        = 0
geometry.is_periodic = 1 1 0

# Keep the one-step fixture runnable on CI nodes with a small host memory
# allocation; the arena grows on demand if the selected physics needs more.
amrex.the_arena_init_size = 33554432

zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

max_step = 1
erf.fixed_dt = 0.1
erf.cfl = 0.5
erf.v = 0
amr.v = 0

erf.plot_file_1 = plt
erf.plot_int_1 = 1
erf.plot_vars_1 = density x_velocity y_velocity z_velocity pressure theta z_phys

erf.init_type = Isentropic
erf.prob_name = "bellforest"
prob.T_0 = 300.0
prob.U_0 = 5.0

erf.use_gravity = true
erf.use_coriolis = false
erf.les_type = "None"
erf.molec_diff_type = "Constant"
erf.dynamic_viscosity = 1.5e-5
erf.alpha_T = 0.0
erf.theta_ref = 300.0
erf.abl_driver_type = "None"

erf.terrain_type = StaticFittedMesh
erf.terrain_smoothing = 2
erf.terrain_file_name_nc = terrain_bell.nc

erf.forest_lai_file = forest_lai_bell.nc
erf.forest_height_file = forest_height_bell.nc
erf.forest_cd = 0.15
erf.forest_tree_type = 2
erf.forest_laimax = 0.6
