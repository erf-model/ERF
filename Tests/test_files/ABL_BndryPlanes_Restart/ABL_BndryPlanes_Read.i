# Reads the BndryFiles series that ABL_BndryPlanes_Restart.i writes across a restart
# (Tests/RunRestartParity.cmake with BNDRY_PLANES_READ_INPUT): the inflow run of
# Exec/CanonicalTests/ABL/inputs_read on the output box [256,768]^2 at the writer's 64 m mesh.
# The harness sets erf.bndry_file to the restarted run's series.
erf.prob_name = "ABL"

max_step = 4

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

erf.substepping_type = None
erf.fixed_dt         = 2.0e-2

erf.v              = 1
amr.v              = 1

amr.max_level       = 0

erf.check_int       = -1
erf.plot_int_1      = -1

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

prob.U_0_Pert_Mag = 0.0
prob.V_0_Pert_Mag = 0.0
prob.W_0_Pert_Mag = 0.0

erf.input_bndry_planes = 1
erf.bndry_file = "BndryFiles"
erf.bndry_input_var_names = temperature density velocity
