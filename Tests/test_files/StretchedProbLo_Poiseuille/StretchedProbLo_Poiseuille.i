# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# PoiseuilleFlow_x on a stretched mesh (uniform ratio) whose bottom is below zero.  The levels
# start at geometry.prob_lo[2] = -1, and z_phys has to follow them: cell centres from -0.9375
# to 0.9375 (erf-model/ERF#4442).
erf.prob_name = "Poiseuille Flow"

max_step = 0

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =  0    0.  -1.
geometry.prob_hi     =  4.   1.   1.
amr.n_cell           = 32    4   16

geometry.is_periodic = 1 1 0

zlo.type = "NoSlipWall"
zhi.type = "NoSlipWall"

# VERTICAL MESH: 16 layers of 0.125 from prob_lo to prob_hi
erf.initial_dz            = 0.125
erf.grid_stretching_ratio = 1.0

# TIME STEP CONTROL
erf.substepping_type   = None
erf.cfl                = 0.5

# DIAGNOSTICS & VERBOSITY
erf.v              = 1
amr.v              = 1

# REFINEMENT / REGRIDDING
amr.max_level       = 0

# CHECKPOINT FILES
erf.check_int       = -1

# PLOTFILES
erf.plot_file_1     = plt
erf.plot_int_1      = 1
erf.plot_vars_1     = x_velocity z_phys

# SOLVER CHOICE
erf.use_gravity            = false

erf.init_type = Uniform

erf.les_type = "None"

erf.molec_diff_type   = "Constant"
erf.dynamic_viscosity = 0.1
erf.alpha_T           = 0.0
erf.alpha_C           = 0.0

erf.abl_driver_type = "PressureGradient"
erf.abl_pressure_grad = -0.2 0. 0.

# PROBLEM PARAMETERS
prob.prob_type = 10

prob.rho_0 = 1.0
prob.T_0 = 300.0
