# Precursor of BndryPlanes_PeriodicSeam: a doubly periodic box with velocity perturbations
# drawn by the cell-indexed hash (erf.deterministic_ic_pert), so the planes it writes vary
# along the faces and are the same on any number of ranks.
erf.prob_name = "ABL"
max_step = 12
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo =    0.    0.     0.
geometry.prob_hi = 1024. 1024.  1024.
amr.n_cell       =   16    16     16
geometry.is_periodic = 1 1 0
zlo.type = "NoSlipWall"
zhi.type = "SlipWall"

erf.timestepping_type = RK3_NoSubstepping
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
erf.deterministic_ic_pert = 1
prob.rho_0 = 1.0
prob.A_0 = 1.0
prob.T_0 = 300.0
prob.U_0 = 10.0
prob.V_0 = 0.0
prob.W_0 = 0.0
prob.U_0_Pert_Mag = 0.5
prob.V_0_Pert_Mag = 0.5
prob.W_0_Pert_Mag = 0.0

erf.output_bndry_planes = 1
erf.bndry_output_planes_interval = 2
erf.bndry_output_start_time = 0.0
erf.bndry_output_planes_file = "BndryFiles"
erf.bndry_output_var_names = temperature velocity density
erf.bndry_output_box_lo = 256. 256.
erf.bndry_output_box_hi = 768. 768.
