# A flat EB covers the physical zlo boundary. Its neutral MOST drag must
# operate with zlo = SlipWall and match the historical zlo = surface_layer run.
erf.prob_name = "ABL"
erf.init_type = "ConstantDensity"
prob.rho_0 = 1.0
prob.T_0 = 300.0
prob.U_0 = 10.0

geometry.prob_lo = 0 0 0
geometry.prob_hi = 32 32 32
amr.n_cell = 16 16 16
geometry.is_periodic = 1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.terrain_type = "EB"
eb2.geometry = "plane"
eb2.plane_point = 0 0 3.25
eb2.plane_normal = 0 0 -1
erf.eb_boundary_type = "SurfaceLayer"
erf.write_eb_surface = false
erf.is_land = 1
erf.most.z0 = 0.1
erf.most.zref = 10.0
erf.most.time_average = true
# Retain history so dry virtual-temperature averaging is exercised.
erf.most.time_window = 10.0
erf.most.surf_temp_flux = 0.0
erf.surface_layer.flux_type = "moeng"

erf.les_type = "Smagorinsky"
erf.Cs = 0.17
erf.molec_diff_type = "None"
erf.use_gravity = false
erf.use_coriolis = false
erf.project_initial_velocity = false
erf.substepping_type = "None"
# Keep the acoustic CFL below one on both levels with no acoustic substepping.
# Spatially varying drag excites modes absent from the uniform-roughness case.
erf.fixed_dt = 0.0005
max_step = 4
erf.check_for_nans = 2
amrex.fpe_trap_invalid = 1
amrex.fpe_trap_zero = 1
amrex.fpe_trap_overflow = 1

amr.max_level = 0
amr.max_grid_size = 8
amr.blocking_factor = 2
amr.n_error_buf = 0
amr.ref_ratio = 2
erf.regrid_int = 1
erf.refinement_indicators = wall
erf.wall.max_level = 1
erf.wall.in_box_lo = 8 8 0
erf.wall.in_box_hi = 24 24 16
erf.wall.start_time = 0.00075
erf.coupling_type = "OneWay"
erf.cf_width = 0
erf.cf_set_width = 0

erf.plot_file_1 = "plt"
erf.plot_int_1 = 4
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta
erf.check_file = "chk"
erf.check_int = 1
erf.v = 1
amr.v = 1
