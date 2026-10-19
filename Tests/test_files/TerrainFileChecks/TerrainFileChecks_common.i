# Shared part of TerrainFileChecks.i (ERF text format) and TerrainFileChecksUSGS.i: checks of
# the terrain text file reader (Source/Utils/ERF_TerrainFileChecks.H) and the warning for steep
# terrain next to an inflow or outflow face.  A compressible flow over a north-south ridge,
# with inflow at xlo and outflow at xhi; make_tables.py writes the tables.
erf.prob_name = "ABL"
erf.init_type = Isentropic

max_step  = 2
stop_time = 1.0e6

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY: 800 m x 400 m, 50 m cells in x and y
geometry.prob_lo     =   0.   0.   0.
geometry.prob_hi     = 800. 400. 800.
amr.n_cell           =  16    8   16
geometry.is_periodic =   0    1    0
amr.max_level        =   0

xlo.type     = "Inflow"
xhi.type     = "Outflow"
xlo.velocity = 5. 0. 0.
xlo.density  = 1.16
xlo.theta    = 300.

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TERRAIN
erf.terrain_type      = StaticFittedMesh
erf.terrain_smoothing = 0

# TIME STEP CONTROL
erf.cfl = 0.5

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = -1
erf.v            = 0
amr.v            = 0

# CHECKPOINT AND PLOT FILES
erf.check_file  = chk
erf.check_int   = -1
erf.plot_file_1 = plt
erf.plot_int_1  = -1
erf.plot_vars_1 = density x_velocity z_phys

# SOLVER CHOICE
erf.use_gravity     = true
erf.les_type        = "None"
erf.pbl_type        = "None"
erf.molec_diff_type = "None"

# PROBLEM PARAMETERS
prob.T_0   = 300.0
prob.U_0   = 5.0
prob.rho_0 = 1.16
