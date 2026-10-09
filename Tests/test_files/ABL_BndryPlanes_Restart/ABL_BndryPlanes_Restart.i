# Boundary-plane output across a restart (Tests/RunRestartParity.cmake with BNDRY_PLANES_DIR).
# The ABL write deck (Exec/CanonicalTests/ABL/inputs_write) on a 16^3 mesh: planes every 2 steps
# from the box [256,768]^2 into BndryFiles. The tests run it straight to step 8 and again through
# a checkpoint, and require the two BndryFiles series to agree plane for plane.
erf.prob_name = "ABL"

max_step = 8

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo =    0.    0.     0.
geometry.prob_hi = 1024. 1024.  1024.
amr.n_cell       =   16    16     16

geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

erf.substepping_type = None
erf.fixed_dt         = 2.0e-2

erf.v              = 1
amr.v              = 1

amr.max_level       = 0

erf.check_file      = chk
erf.check_int       = -1

erf.plot_file_1     = plt
erf.plot_int_1      = -1
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

prob.U_0_Pert_Mag = 0.08
prob.V_0_Pert_Mag = 0.08
prob.W_0_Pert_Mag = 0.0

erf.output_bndry_planes = 1
erf.bndry_output_planes_interval = 2
erf.bndry_output_start_time = 0.0
erf.bndry_output_planes_file = "BndryFiles"
erf.bndry_output_var_names = temperature velocity density

erf.bndry_output_box_lo = 256. 256.
erf.bndry_output_box_hi = 768. 768.
