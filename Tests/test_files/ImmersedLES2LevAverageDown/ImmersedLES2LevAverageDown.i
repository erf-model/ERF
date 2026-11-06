# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Immersed-forcing terrain (a plane at 50 m, the mid-plane of a 20 m coarse cell) with Smagorinsky
# LES, anelastic, two levels (20 m over 10 m, the fine level over the whole width up to 100 m), and a
# stable theta gradient (0.02 K/m), so the solid fine cells hold another theta than the fluid ones.
# With erf.if_average_down = fluid_mass every covered coarse cell takes the fluid-mass-weighted
# average of the fine cells over it; the plain average (the default) mixes the solid fine cell
# (40-50 m) into the coarse cell 40-60 m, about 0.1 K off here.
#
# Tests/CTestList.cmake runs it with both settings and checks the covered coarse cells against the
# fine level (Tests/ImmersedAverageDownCheck.cpp): the fluid_mass run must pass, the plain run must
# not.
#
erf.prob_name = "ABL"

max_step = 5

amrex.fpe_trap_invalid  = 1
amrex.fpe_trap_zero     = 1
amrex.fpe_trap_overflow = 1

geometry.prob_lo     = 0.0 0.0 0.0
geometry.prob_hi     = 320.0 160.0 400.0
amr.n_cell           = 16 8 20
amr.max_grid_size    = 8        # several boxes per level, so the fine sums cross ranks
amr.blocking_factor  = 2
geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"   # inside the solid; the immersed forcing carries the wall
zhi.type = "SlipWall"

erf.terrain_type      = ImmersedForcing
erf.terrain_file_name = "terrain.txt"
erf.if_use_most       = true
erf.if_z0             = 0.1
erf.if_implicit_drag  = true
erf.if_average_down   = fluid_mass
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
erf.plot_int_1  = 5
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
erf.box1.in_box_lo_indices_crse = 0 0 0
erf.box1.in_box_hi_indices_crse = 15 7 4
