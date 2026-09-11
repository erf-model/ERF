# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Two-level BTF terrain-following mesh over a radial Witch of Agnesi hill, used by the
# Terrain2Lev_BTF_ZSplit parity test (Tests/RunTerrainZSplitParity.cmake).
#
# The test runs this deck twice to step 0: once with whole-height fine grids
# (amr.max_grid_size_z = 1024 1024) and once with the fine grids split at k = 16 and 32
# (amr.max_grid_size_z = 1024 16), and requires the plotfiles to agree bit for bit.
# The fine level covers the hill top, where the fine terrain differs most from the
# terrain interpolated from the coarse level.
#
# Only the mesh is plotted: on a fine level the hydrostatic base state is still
# integrated box by box (init_bcs_and_base_state), so dens_hse and pres_hse there still
# depend on the z split.

erf.prob_name = "Flow over Witch of Agnesi hill"

erf.init_type = Isentropic

max_step = 0

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =   0.   0.   0.
geometry.prob_hi     = 640. 640. 400.
amr.n_cell           =  32   32   40

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 0.1

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = -1
erf.v            = 0
amr.v            = 0

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 1
erf.plot_vars_1 = z_phys detJ

# SOLVER CHOICE
erf.use_gravity     = true
erf.molec_diff_type = "None"
erf.les_type        = "None"

# TERRAIN GRID TYPE
erf.terrain_type      = StaticFittedMesh
erf.terrain_smoothing = 0                       # BTF

# MULTILEVEL
amr.max_level      = 1
amr.ref_ratio_vect = 2 2 1

erf.refinement_indicators = box1
erf.box1.max_level = 1
erf.box1.in_box_lo = 160. 160.
erf.box1.in_box_hi = 480. 480.

amr.max_grid_size_x = 1024 16
amr.max_grid_size_y = 1024 16
amr.max_grid_size_z = 1024 1024

# PROBLEM PARAMETERS
prob.T_0   = 300.0
prob.U_0   = 0.0
prob.rho_0 = 1.16

prob.custom_terrain_type = "WoA"
prob.dir                 = 2                    # radial hill
prob.hmax                = 100.0
prob.L                   = 100.0
