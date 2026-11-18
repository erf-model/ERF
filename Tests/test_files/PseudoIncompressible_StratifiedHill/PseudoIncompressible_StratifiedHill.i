# Pseudo-incompressible flow over the 100 m hill of AnelasticScalarAdvection/inputs_hill, in a
# stable sounding (theta = 300 + 0.03 z K), on two levels coupled two-way.  The density and
# (rho theta) evolve, are refluxed and are averaged down, so their composite integrals must
# stay the same (Tests/RunConservation.cmake).

erf.prob_name = "Scalar Advection/Diffusion"
stop_time = 1.0e6
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024
geometry.prob_lo     = -1280. -1280.    0.
geometry.prob_hi     =  1280.  1280.  791.6  # 12 levels from 20 m stretched by 1.2
amr.n_cell           =    32     32   12     # dx = dy = 80 m
geometry.is_periodic = 1 1 0
erf.terrain_type          = StaticFittedMesh
erf.terrain_smoothing     = 0
erf.initial_dz            = 20.0
erf.grid_stretching_ratio = 1.2
prob.custom_terrain_type  = "WoA"
prob.dir                  = 2
prob.hmax                 = 100.0
prob.L                    = 300.0
erf.terrain_poisson_solver = mlmg
zlo.type = "SlipWall"
zhi.type = "SlipWall"
erf.fixed_dt = 2.0        # advective CFL about 0.3
amr.v = 0
erf.plot_file_1 = plt
erf.plot_vars_1 = density scalar theta x_velocity y_velocity z_velocity z_phys
erf.dycore_horiz_adv_type    = "Upwind_3rd"
erf.dycore_vert_adv_type     = "Upwind_3rd"
erf.dryscal_horiz_adv_type   = "Upwind_3rd"
erf.dryscal_vert_adv_type    = "Upwind_3rd"
erf.les_type        = "None"
erf.molec_diff_type = "None"
prob.rho_0 = 1.0
prob.T_0   = 300.0
prob.A_0   = 1.0
prob.U_0   = 10.0
prob.V_0   = 3.0
prob.rad_0 = 100.0        # a ball far larger than the domain: scalar = A_0 everywhere

# ---- Overrides for the pseudo-incompressible conservation test
erf.compressibility_type = PseudoIncompressible
erf.use_gravity          = true
erf.init_type            = input_sounding
erf.input_sounding_file  = sounding_stable
prob.test_flow_amp       = 2
max_step                 = 10
erf.plot_int_1           = -1
erf.check_int            = -1
erf.v                    = 1
erf.sum_interval         = 1
amr.max_level             = 1
amr.ref_ratio_vect        = 2 2 1
erf.refinement_indicators = box1
erf.box1.max_level        = 1
erf.box1.in_box_lo_indices_crse = 8 12 0
erf.box1.in_box_hi_indices_crse = 19 23 11
