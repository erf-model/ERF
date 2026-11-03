# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# The flat.i ABL over an immersed plane on the cell face z = 0, two solid rows below it; otherwise as flat.i
# (see its header). Flat neutral ABL with the surface layer at the bottom (k-eqn RANS, anelastic MidPoint, implicit
# vertical diffusion, FFT Poisson solve): the reference for plane.i, the same ABL over an immersed
# plane that lies on a cell face, with the fraction-stress wall law (erf-model/ERF#4134). On a
# face-aligned surface the fraction-stress law is the surface-layer law of flat ground at dz/2, the
# wall face is closed to diffusion like zlo, and the RANS wall distance is measured from the plane,
# so the two runs agree cell by cell (shifted by the two solid rows: wind within 0.02 %, TKE within
# 0.2 % over the lowest 500 m after 100 steps); the legacy law leaves
# the fluid cell on the solid unforced and measures the wall distance from the domain bottom, and
# must fail the same check. Tests/RunImmersedFlatEquivalence.cmake runs both decks and compares the
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
erf.anelastic      = 1
erf.anelastic_type = MidPoint
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
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta KE Kmv walldist

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
