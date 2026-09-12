# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Anelastic neutral flow over a radial Witch of Agnesi hill on a terrain-fitted
# mesh under a MOST surface layer, with more cells in z than amr.max_grid_size so
# that every column is several boxes. That box split used to break (1) the initial
# projection, which read the momenta one ghost face beyond each box face inside
# the domain before those faces were filled, and (2) the planar surface-layer
# arrays (u*, t*, ...), which hold one 2D box per 3D box and were FillBoundary'ed
# across the uncomputed duplicate copies. The domain is also split into several
# columns, so the fill of the surface copies crosses box boundaries in the plane.
erf.prob_name = "ABL"

max_step = 10

erf.anelastic = 1

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -1280. -1280.   0.
geometry.prob_hi     =  1280.  1280. 800.
amr.n_cell           =    32     32   40     # dx = dy = 80 m, dz = 20 m
amr.max_grid_size    = 16                    # 4 columns of 3 boxes, split in z at k = 16 and 32
amr.blocking_factor  = 8

geometry.is_periodic = 1 1 0

# TERRAIN: axisymmetric Witch of Agnesi hill centred in the domain
erf.terrain_type         = StaticFittedMesh
erf.terrain_smoothing    = 0
prob.custom_terrain_type = "WoA"
prob.dir                 = 2        # radial
prob.hmax                = 100.0
prob.L                   = 300.0

# BOUNDARY CONDITIONS
zlo.type      = "surface_layer"
erf.most.z0   = 0.1
erf.most.zref = 10.0

zhi.type      = "SlipWall"

# INITIALIZATION
erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"

# TIME STEP CONTROL
erf.fixed_dt = 1.5

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.v            = 1
amr.v            = 1

# REFINEMENT / REGRIDDING
amr.max_level = 0

# CHECKPOINT FILES
erf.check_file = chk
erf.check_int  = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 10
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta

# SOLVER CHOICE
erf.use_gravity  = true
erf.use_coriolis = false

erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.1

# PROBLEM PARAMETERS: no perturbations
prob.T_0_Pert_Mag = 0.0
prob.U_0_Pert_Mag = 0.0
prob.V_0_Pert_Mag = 0.0
prob.W_0_Pert_Mag = 0.0
