# InflowThetaDensity with xlo.type = inflow_outflow: on an ext_dir_upwind face a primitive
# theta from the dirichlet_file must still be multiplied by the ghost density. Assigned to
# RhoTheta directly it is an effective theta of 300 / 0.95 = 316 K. The test requires theta
# to stay 300 K.
erf.prob_name = "ABL"
max_step = 10
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_extent =   500     500    500
amr.n_cell           =    32      32     32
geometry.is_periodic = 0 1 0

xlo.type = "inflow_outflow"
xlo.density = 0.95
xlo.dirichlet_file = "inflow_file"
xhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.fixed_dt = 0.1
erf.fixed_mri_dt_ratio = 6
erf.sum_interval = 1
erf.v   = 1
amr.v   = 1
amr.max_level = 0

erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 10
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta

erf.use_gravity     = false
erf.molec_diff_type = "None"
erf.les_type        = "None"

erf.init_type = Uniform
prob.rho_0 = 1.0
prob.T_0   = 300.0
prob.U_0   = 10.0
prob.V_0   = 0.0
prob.W_0   = 0.0
