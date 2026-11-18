# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# The density current of Tests/test_files/DensityCurrent at dx = dz = 200 m, as the field that
# the refinement tests of Tests/check_tagging_conditions.py tag on: the cold bubble gives theta,
# and density varies with both the bubble and the height. The refinement indicator itself is
# given on the command line by each test.
erf.prob_name = "Density Current"

erf.init_type = Isentropic

max_step = 0
stop_time = 900.0

erf.buoyancy_type = 1

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -12800.   0.    0.
geometry.prob_hi     =  12800. 100. 6400.
amr.n_cell           =  128      4    32

geometry.is_periodic = 0 1 0

xlo.type = "Symmetry"
xhi.type = "Outflow"

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt       = 1.0
erf.fixed_fast_dt  = 0.25

# DIAGNOSTICS & VERBOSITY
erf.v              = 1
amr.v              = 1

# REFINEMENT / REGRIDDING
# With no error buffer, a coarse blocking factor of 1 and a grid efficiency of 1, the level 1
#    grids cover exactly the tagged cells
amr.max_level       = 0
amr.ref_ratio       = 2
amr.blocking_factor = 2
amr.max_grid_size   = 32
amr.n_error_buf     = 0
amr.grid_eff        = 1.0

# CHECKPOINT FILES
erf.check_int       = -1

# PLOTFILES
erf.plot_file_1     = plt
erf.plot_int_1      = 1
erf.plot_vars_1     = density theta

# SOLVER CHOICE
erf.use_gravity = true

erf.les_type         = "None"

# The level 1 grids follow the tags, so they may be split in z
erf.vert_implicit_fac = 0 0 0

erf.molec_diff_type   = "ConstantAlpha"
erf.dynamic_viscosity = 75.0
erf.alpha_T           = 75.0
erf.alpha_C           = 0.0

erf.c_p = 1004.0

# PROBLEM PARAMETERS
prob.T_0 = 300.0
