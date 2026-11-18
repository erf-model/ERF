# BndryPlanes_PeriodicSeam: planes ingested on the x faces of a y-periodic domain. The planes
# vary along y (the writer's perturbations), so the ghost cells beside the x faces at the y seam
# must hold the plane at the periodic image (j = ny - 1 for j = -1), as the interior cells across
# the seam do; development's near-side fill changed y_velocity by 2.7 % of its norm after
# 10 steps (6.6 % on the canonical ABL inputs_write / inputs_read pair). The velocity components
# nodal along a face are read between the two plane cells either side of the node, which the plane
# holds at cell centres (1.7 % in y_velocity here). The reader's last
# plotfile is compared with the gold, which records the new fill, not a physical answer: the
# reader's 512 m window is not a period of the writer's 1024 m field, so neither fill is "the"
# value at the seam; the BndryPlaneIndex gtests carry the check of the index arithmetic.
erf.prob_name = "ABL"
max_step = 10
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo =  256.  256.     0.
geometry.prob_hi =  768.  768.  1024.
amr.n_cell       =    8     8     16
geometry.is_periodic = 0 1 0
xlo.type = "Inflow"
xhi.type = "Outflow"
zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

erf.timestepping_type = RK3_NoSubstepping
erf.fixed_dt         = 2.0e-2
erf.v              = 1
amr.v              = 1
amr.max_level       = 0
erf.check_int       = -1
erf.plot_file_1     = plt
erf.plot_int_1      = 10
erf.plot_vars_1     = density x_velocity y_velocity z_velocity theta

erf.use_gravity = false
erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.1

erf.init_type = "uniform"
prob.rho_0 = 1.0
prob.A_0 = 1.0
prob.T_0 = 300.0
prob.U_0 = 10.0
prob.V_0 = 0.0
prob.W_0 = 0.0

erf.input_bndry_planes = 1
erf.bndry_file = "BndryFiles"
erf.bndry_input_var_names = temperature density velocity
