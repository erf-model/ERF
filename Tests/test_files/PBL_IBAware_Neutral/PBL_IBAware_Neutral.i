# Decomposition parity of YSUNew with PBL-height smoothing over a neutral
# surface (erf-model/ERF#4479), based on inputs_base of
# Exec/RegTests/ImmersedForcingTest/PBL_IBAware (one 32 x 32 box, YSUNew,
# 10 steps, 2D output added), run tiled and untiled by
# Tests/RunTilingParity.cmake and on one box and on 16 by
# Tests/RunBoxParity.cmake. The ground is neutral (300 K under a 300 K
# sounding), so the surface heat flux is round-off whose sign depends on the
# decomposition. The bulk Richardson number reaches Ribcr in no column, so
# every column takes the top of the scan, 0.9 z(khi) = 139.5 m (#4458), and
# the diffusivity through that depth is the one the surface stability flag
# selects. pblh is the same in every column, so the tiling test checks only
# that Kmv and Khv vary; restoring the old flag rule (sflux > 0) still makes
# both decomposition tests fail.
stop_time = 100000.0
amrex.fpe_trap_invalid = 1

geometry.prob_extent =  320     320    160
amr.n_cell           =  32      32     16
amr.max_grid_size    =  32
geometry.is_periodic = 1 1 0
zlo.type = "surface_layer"
erf.most.z0 = 0.1
erf.most.pot_surf_temp = 300.0
erf.most.pblh_calc = "YSU"
zhi.type = "SlipWall"
erf.fixed_dt           = 0.5
erf.substepping_cfl    = 0.5
erf.sum_interval   = -1
erf.v              = 0
amr.v              = 0
amr.max_level      = 0
erf.check_int       = -1
erf.plot_file_1     = plt
erf.plot_int_1      = 10
erf.plot_vars_1     = theta x_velocity terrain_IB_mask Kmv Khv Lturb
erf.plot2d_file_1   = plt2d
erf.plot2d_int_1    = 10
erf.plot2d_vars_1   = pblh u_star t_star
erf.use_gravity = true
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.buildings_type = ImmersedForcing
eb2.geometry = box
eb2.box_lo = 140.0 140.0 -10.0
eb2.box_hi = 180.0 180.0 40.0
erf.immersed_forcing_substep = true
eb2.small_volfrac = 0.005
erf.if_use_most = true
erf.if_z0 = 0.01
erf.use_coriolis = false
erf.pbl_type = YSUNew
erf.pbl_ib_aware = true
erf.pbl_ib_z0 = 0.01
erf.enable_pblh_smoothing = true
erf.pblh_smoothing_passes = 1
max_step = 10
