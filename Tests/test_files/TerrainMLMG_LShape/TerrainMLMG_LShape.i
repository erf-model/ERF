# Two-level anelastic run over a radial hill whose refined level is an L-shaped
# union of two boxes.  The FFT-preconditioned GMRES projection needs a rectangular
# union of boxes on every level and stops here; the multigrid terrain solver
# (erf.terrain_poisson_solver = mlmg or gmres_mlmg) runs.  Coarsened from
# Exec/CanonicalTests/Canonical_RANS/Neutral_Hill_3D.
erf.prob_name = "ABL"

erf.anelastic = 1
erf.terrain_poisson_solver = mlmg
erf.mg_v = 1            # prints the divergence before and after every projection

stop_time = 100.0
amrex.fpe_trap_invalid = 0

geometry.prob_lo     = -2560. -2560.   0.
geometry.prob_hi     =  2560.  2560. 800.
amr.n_cell           =    32     32   16     # dx = dy = 160 m, dz = 50 m
amr.max_grid_size    = 16
amr.blocking_factor  = 4
geometry.is_periodic = 1 1 0

# Refined level: two boxes whose union is not a rectangle.  The regridder must keep
# that shape rather than fill it in to its bounding box, so no buffer cells and a
# grid efficiency of one.
amr.max_level = 1
amr.grid_eff = 1.0
amr.n_error_buf = 0
erf.coupling_type = OneWay
erf.refinement_indicators = box_south box_northeast
erf.box_south.max_level     = 1
erf.box_south.in_box_lo     = -1280. -1280.
erf.box_south.in_box_hi     =  1280.     0.
erf.box_northeast.max_level = 1
erf.box_northeast.in_box_lo =     0.     0.
erf.box_northeast.in_box_hi =  1280.  1280.

erf.terrain_type         = StaticFittedMesh
erf.terrain_smoothing    = 0
erf.wall_dist_type       = terrain_height
prob.custom_terrain_type = "WoA"
prob.dir                 = 2        # radial
prob.hmax                = 100.0    # hill height [m]
prob.L                   = 500.0    # half-width at half height [m]

zlo.type      = "surface_layer"
erf.most.z0   = 0.1
zhi.type       = "SlipWall"
zhi.theta_grad = 0.03

erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.sounding_type       = "Ideal"

erf.fixed_dt = 1.5

erf.sum_interval = 10
erf.v            = 1
amr.v            = 0

erf.check_file = chk
erf.check_int  = -1
erf.plot_file_1 = plt
erf.plot_int_1  = -1
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta KE z_phys

erf.dycore_horiz_adv_type  = "Upwind_5th"
erf.dycore_vert_adv_type   = "Upwind_3rd"
erf.dryscal_horiz_adv_type = "Upwind_5th"
erf.dryscal_vert_adv_type  = "Upwind_3rd"

erf.use_gravity  = true
erf.use_coriolis = true
erf.coriolis_3d  = false
erf.latitude     = 90.0
erf.rotational_time_period = 125663.7061435917   # f = 1e-4 1/s
erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind    = 10.0 0.0 0.0

erf.rayleigh_damp_W   = true
erf.rayleigh_dampcoef = 0.2
erf.rayleigh_zdamp    = 250.

erf.molec_diff_type = "None"
erf.les_type            = "None"
erf.rans_type           = "kEqn"
erf.dirichlet_k         = true
erf.init_tke_from_ustar = true
erf.theta_ref           = 300.0

prob.pert_deltaU  = 0.0
prob.pert_deltaV  = 0.0
prob.T_0_Pert_Mag = 0.0
prob.U_0_Pert_Mag = 0.0
prob.V_0_Pert_Mag = 0.0
prob.W_0_Pert_Mag = 0.0
