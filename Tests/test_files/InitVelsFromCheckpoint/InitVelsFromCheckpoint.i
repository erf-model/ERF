# Small ABL deck for the erf.init_vels_from_checkpoint test (issue 5376).
# RunInitVelsFromCheckpoint.cmake runs it to chk00005 on 16^3 boxes, then starts it again
# on other boxes with erf.init_vels_from_checkpoint = chk00005.
erf.prob_name = "ABL"

max_step = 5

amrex.fpe_trap_invalid = 1

geometry.prob_extent = 1024 1024 512
amr.n_cell           =   32   32  16
amr.max_grid_size    = 16

geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

erf.fixed_dt = 0.1

erf.v = 1
amr.v = 1

amr.max_level = 0

erf.check_file = chk
erf.check_int  = 5

erf.plot_file_1 = plt
erf.plot_int_1  = 5
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta

erf.use_gravity = false

erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.1

erf.init_type = "uniform"

prob.rho_0 = 1.0
prob.A_0 = 1.0

prob.U_0 = 10.0
prob.V_0 = 0.0
prob.W_0 = 0.0
prob.T_0 = 300.0

prob.U_0_Pert_Mag = 0.08
prob.V_0_Pert_Mag = 0.08
prob.W_0_Pert_Mag = 0.0
