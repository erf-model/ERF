# Tiling parity of the immersed-boundary-aware MRF and YSUNew schemes
# (erf.pbl_ib_aware): a 40 m cube by immersed forcing under a surface layer,
# run tiled and untiled by Tests/RunTilingParity.cmake. The scheme is set
# from the CTest entry. The sounding caps a 60 m mixed layer with an
# inversion so the diagnosed boundary-layer height varies from column to
# column, and is lower above the 40 m roof than over open ground.
# The MRF and YSUNew boundary-layer schemes on an immersed building.
#
# A 40 m cube (an exact box, so every cell is fluid or solid) on a 320 m
# periodic domain at 10 m in a 3 m/s westerly, neutral at 300 K with a MOST
# ground, the PBL scheme supplying the vertical diffusivity. The schemes
# measure height from the terrain surface, so over the cube the wall
# distance and the boundary-layer depth are wrong by the building height and
# the diffusivity profile continues inside the solid. With
# erf.pbl_ib_aware the column's surface is the first fluid cell above the
# solid: heights and the depth are measured from the roof, the diffusivity
# is zero inside the solid, and the surface scales over the cube come from a
# neutral log law at the roof with erf.pbl_ib_z0.
# Variants: inputs_mrf_off/on, inputs_ysunew_off/on, and inputs_flat_off/on
# (no building: the switch must change nothing).
stop_time = 100000.0

geometry.prob_extent =  320     320    160
amr.n_cell           =  32      32     16
amr.max_grid_size    =  16
geometry.is_periodic = 1 1 0
zlo.type = "surface_layer"
erf.most.z0 = 0.1
erf.most.surf_temp_flux = 0.1     # heated surface, so the PBL height varies from column to column
erf.most.pblh_calc = "MRF"
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
erf.plot2d_vars_1   = pblh u_star
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
erf.pbl_ib_z0 = 0.01
erf.pbl_ib_aware = true
erf.pbl_type = "MRF"
max_step = 10
