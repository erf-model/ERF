# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Immersed-forcing terrain (a plane at 45 m) with Smagorinsky LES, anelastic, two levels: a refined
# patch (20 m over 10 m) whose lateral coarse/fine faces cross the surface. Theta is uniform and no
# heat enters, so no net heat may appear (erf-model/ERF#4920). With the implicit projection
# (erf.if_implicit_projection = true, erf-model/ERF#4926) the solid is not frozen after the
# projection and the momenta stay divergence free in the partial cells at the patch edges. With
# the explicit freeze (erf.if_implicit_projection = false) those partial cells keep a divergence
# of their own, which advects the uniform theta. Measured on one build (Release, 2 ranks, 20
# steps): mean theta drift 2e-9 K implicit, 0.04 K with the freeze (1 K in the half-solid fine
# layer at the inflow edge); one level alone 1e-7 K either way.
#
# Tests/CTestList.cmake runs it with both settings and checks the mean theta of the domain, the
# fine level where it exists (Tests/ImmersedCFHeatBalanceCheck.cpp): the implicit run must hold
# it within 1e-6 K, the freeze run must not. A third leg adds the log-law coarse/fine fill
# (erf.cf_loglaw_fill) and a mismatch threshold the surface crosses (erf.if_cf_mismatch_warning =
# 0.4: the fine cell 40-50 m is half solid, a spread of 0.5) and must hold the mean and print the
# warning. The implicit projection needs an FFT build.
#
erf.prob_name = "ABL"

max_step = 20

amrex.fpe_trap_invalid  = 1
amrex.fpe_trap_zero     = 1
amrex.fpe_trap_overflow = 1

geometry.prob_lo     = 0.0 0.0 0.0
geometry.prob_hi     = 320.0 160.0 400.0
amr.n_cell           = 16 8 20
amr.max_grid_size    = 256
amr.blocking_factor  = 2
geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"   # inside the solid; the immersed forcing carries the wall
zhi.type = "SlipWall"

erf.terrain_type      = ImmersedForcing
erf.terrain_file_name = "terrain.txt"
erf.if_use_most       = true
erf.if_z0             = 0.1
erf.if_implicit_drag  = true
erf.if_implicit_projection = true
eb2.cover_multiple_cuts = 1

erf.compressibility_type = Anelastic
erf.timestepping_type = MidPoint
erf.cfl            = 0.5

erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.sounding_type       = "Ideal"

erf.v = 1
amr.v = 0

erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 20
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta terrain_IB_mask

erf.dycore_horiz_adv_type  = "Upwind_3rd"
erf.dycore_vert_adv_type   = "Upwind_3rd"
erf.dryscal_horiz_adv_type = "Upwind_3rd"
erf.dryscal_vert_adv_type  = "Upwind_3rd"

erf.use_gravity     = true
erf.use_coriolis    = true
erf.coriolis_3d     = false
erf.latitude        = 90.0
erf.rotational_time_period = 125663.7061435917
erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind    = 10.0 0.0 0.0

erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.17

prob.pert_deltaU  = 0.0
prob.pert_deltaV  = 0.0
prob.T_0_Pert_Mag = 0.0

amr.max_level     = 1
amr.ref_ratio     = 2 2 2
erf.coupling_type = "TwoWay"
erf.refinement_indicators = box1
erf.box1.max_level = 1
erf.box1.in_box_lo_indices_crse = 4 0 0
erf.box1.in_box_hi_indices_crse = 11 7 4
