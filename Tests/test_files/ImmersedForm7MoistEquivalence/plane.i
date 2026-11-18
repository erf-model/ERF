# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# The flat.i moist ABL over an immersed plane on the cell face z = 0, two dry solid rows (0 g/kg)
# below it; otherwise as flat.i (see its header). With the fraction-stress wall law the wall face
# is closed to the moisture species and they are held in the solid, so the qv profile over the
# plane matches flat ground; a wall face open to moisture dries the wall cell by tens of percent
# within the 100 steps. Tests/RunImmersedFlatEquivalence.cmake runs both decks and compares the
# horizontally averaged profiles with Tests/ImmersedFlatEquivalenceCheck.cpp.
#
erf.prob_name = "ABL"
stop_time = 1.0e6
max_step  = 100
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo     = 0.0 0.0 -40.0
geometry.prob_hi     = 160.0 160.0 1000.0
amr.n_cell           = 8 8 52
amr.max_grid_size    = 256
amr.blocking_factor  = 2
amr.max_level        = 0
geometry.is_periodic = 1 1 0

zlo.type    = "NoSlipWall"   # inside the solid; the immersed forcing carries the wall
# (erf.dirichlet_k needs the surface layer; the wall TKE comes from erf.if_wall_tke below)
zhi.type       = "SlipWall"
zhi.theta_grad = 0.0

# anelastic, implicit vertical diffusion (solved in the first MidPoint stage)
erf.compressibility_type = Anelastic
erf.timestepping_type = MidPoint
erf.vert_implicit  = true
erf.use_fft        = true    # FFT Poisson solve (flat, periodic in x and y)
erf.fixed_dt       = 1.0

erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding_plane"
erf.sounding_type       = "Ideal"

erf.sum_interval    = 60
erf.v               = 1
amr.v               = 0

erf.check_file  = chk
erf.check_int   = -1
erf.plot_file_1 = plt_plane
erf.plot_int_1  = 100
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta KE Kmv walldist qv

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

# top sponge toward the sounding (10 m/s), as the kynema-sgf deck; dampcoef * dt = 0.05
erf.rayleigh_damp_U   = true
erf.rayleigh_damp_V   = true
erf.rayleigh_damp_W   = true
erf.rayleigh_dampcoef = 0.05
erf.rayleigh_zdamp    = 200.0

erf.moisture_model  = "MoistNoCondensation"   # qv and qc transported, no condensation
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.rans_type       = "kEqn"
erf.wall_dist_type  = terrain_height   # legacy law: the distance from the domain bottom (no effect with fraction_stress)

erf.terrain_type      = ImmersedForcing
erf.terrain_file_name = "terrain.txt"
erf.if_use_most       = true
erf.if_z0             = 0.1
erf.if_implicit_drag  = true
erf.if_wall_form      = fraction_stress
erf.if_wall_tke       = true
eb2.cover_multiple_cuts = 1
erf.theta_ref       = 300.0

# same initial k in every case: 0.84 m2/s2, linear decay to zero at 900 m above z = 0
# the initial TKE profile is measured from geometry.prob_lo, 40 m below the ground: KE_0 and the
# decay height are scaled so that it matches flat.i's 0.84 (1 - z/900) above the ground
prob.KE_0            = 0.877333333
prob.KE_decay_height = 940.0
prob.KE_decay_order  = 1
prob.pert_deltaU  = 0.0
prob.pert_deltaV  = 0.0
prob.T_0_Pert_Mag = 0.0
prob.U_0_Pert_Mag = 0.0
prob.V_0_Pert_Mag = 0.0
prob.W_0_Pert_Mag = 0.0
