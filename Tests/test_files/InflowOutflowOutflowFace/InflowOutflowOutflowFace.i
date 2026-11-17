# Flow leaves through an inflow_outflow face (xhi; xlo is an Inflow face at 12 m/s into a
# 10 m/s initial flow, so the state at xhi changes during the run).  Where the flow leaves,
# the ghost cells of the velocities and, with xhi.density given, of the conserved variables
# must be filled by foextrap (#5136).  They used to be left unfilled: a fresh run carried
# their initial values along, and a restart, which sets them to bogus_large_value, read
# those.  The test requires a run restarted at step 15 to match a straight run exactly.

# ------------------  INPUTS TO MAIN PROGRAM  -------------------
erf.prob_name = "ABL"
max_step = 30

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =   500     500    500
amr.n_cell           =    32      32     32

geometry.is_periodic = 0 1 0

xlo.type     = "Inflow"
xlo.velocity = 12. 0. 0.
xlo.density  = 1.0
xlo.theta    = 300.

xhi.type     = "inflow_outflow"
xhi.velocity = 10. 0. 0.
xhi.density  = 1.0
xhi.theta    = 300.

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt           = 0.1
erf.fixed_mri_dt_ratio = 6

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.v   = 1
amr.v   = 1

# REFINEMENT / REGRIDDING
amr.max_level = 0

# PLOTFILES
erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = 30
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta

# SOLVER CHOICE
erf.use_gravity     = false
erf.molec_diff_type = "None"
erf.les_type        = "None"

erf.init_type = Uniform

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.T_0   = 300.0
prob.U_0   = 10.0
prob.V_0   = 0.0
prob.W_0   = 0.0
