# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Immersed-forcing building (a box 60 x 60 m in plan and 70 m tall, eb2.geometry = box) with
# Smagorinsky LES, anelastic, one level split into four boxes of unequal size (12 + 4 cells in x,
# 4 + 4 in y) whose seams cross the building. The box edges sit on cell mid-planes in x, y and z,
# so its walls and roof are rows of half-solid cells. Theta is uniform and no heat enters, so no
# net heat may appear
# (erf-model/ERF#4887). With the implicit projection (erf.if_implicit_projection = true) the linear
# drag of the building kernels acts inside the projection and the solid is not frozen after it,
# so the momenta stay divergence free in the partial cells at the walls and roof. With the explicit
# treatment (erf.if_implicit_projection = false) the faces between two solid cells are zeroed after
# the projection and the partial cells next to them keep a divergence of their own, which advects
# the uniform theta.
#
# Tests/CTestList.cmake runs it with both settings and checks the fluid-mass-weighted mean theta of
# the domain (Tests/ImmersedCFHeatBalanceCheck.cpp): the implicit run must hold it within 1e-6 K
# (measured on one build: 4e-10 K, and with the snap 2e-9 K), the explicit run must not (measured
# -0.012 K, with theta between 263 and 358 K in the partial cells). A third leg runs the implicit
# projection with erf.if_snap_partial_cells. The implicit projection needs an FFT build.
#
erf.prob_name = "ABL"

max_step = 20

amrex.fpe_trap_invalid  = 1
amrex.fpe_trap_zero     = 1
amrex.fpe_trap_overflow = 1

geometry.prob_lo     = 0.0 0.0 0.0
geometry.prob_hi     = 320.0 160.0 400.0
amr.n_cell           = 16 8 20
amr.max_grid_size_x  = 12
amr.max_grid_size_y  = 4
amr.max_grid_size_z  = 20
amr.blocking_factor  = 2
geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

erf.buildings_type   = ImmersedForcing
eb2.geometry         = box
eb2.box_lo           = 130.0 50.0 -20.0
eb2.box_hi           = 190.0 110.0 70.0
erf.if_use_most      = true
erf.if_z0            = 0.1
erf.if_implicit_drag = true
erf.if_implicit_projection = true

erf.compressibility_type = Anelastic
erf.timestepping_type = MidPoint
erf.use_fft        = true
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

amr.max_level = 0
