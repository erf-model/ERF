# The immersed-boundary-aware MRF scheme (erf.pbl_ib_aware) with the PBL
# height smoothing (erf.enable_pblh_smoothing) on: the same 40 m cube, heated
# ground and capped sounding as PBL_IBAware_Tiling, compared with a gold
# plotfile. The smoothing stencil mixes columns whose surfaces differ by the
# building height, so it must run on the absolute height and keep every
# column at or above its own floor; Lturb (the stored height) in the 3D
# plotfile and the flow above the roof carry the result.
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
erf.enable_pblh_smoothing = true
max_step = 10
