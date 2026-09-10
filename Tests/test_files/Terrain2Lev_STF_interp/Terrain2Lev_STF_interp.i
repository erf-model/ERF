# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Two-level STF terrain-following mesh, "interpolate" refinement mode.
#
# The level-1 grids deliberately touch the non-periodic xlo boundary: that is the
# configuration in which the nodes outside the domain matter, since they are the ones
# make_terrain_fitted_coords poisons on entry and the fine-level path has to refill.
# They also stop below the domain top, so the fine mesh has a coarse/fine interface
# both laterally (at x = 800) and above it.
#
# z_phys and detJ are in the plotfile so that the comparison pins the mesh itself, not
# just the state carried on it.

erf.prob_name = "Flow over Witch of Agnesi hill"

erf.init_type = Isentropic

max_step = 10

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -3200.   0.     0.
geometry.prob_hi     =  3200. 400.  3200.
amr.n_cell           =     64    4     32

geometry.is_periodic = 0 1 0

xlo.type = "Inflow"
xhi.type = "Outflow"
xlo.velocity = 10. 0. 0.
xlo.density  = 1.16
xlo.theta    = 300.
xlo.scalar   = 0.

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt           = 0.5
erf.fixed_mri_dt_ratio = 6

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = -1
erf.v              = 0
amr.v              = 0

# CHECKPOINT FILES
erf.check_file      = chk
erf.check_int       = -1

# PLOTFILES
erf.plot_file_1     = plt
erf.plot_int_1      = 10
erf.plot_vars_1     = density x_velocity y_velocity z_velocity pressure theta z_phys detJ

# SOLVER CHOICE
erf.use_gravity  = true
erf.use_coriolis = false

erf.molec_diff_type = "None"
erf.les_type        = "None"

# TERRAIN GRID TYPE
erf.terrain_type      = StaticFittedMesh
erf.terrain_smoothing = 1                       # STF

# MULTILEVEL
amr.max_level      = 1
amr.ref_ratio_vect = 2 1 2

# Regrid so that remake_zphys() rebuilds the fine mesh at least once
erf.regrid_int = 5

erf.refinement_indicators = box1
erf.box1.max_level = 1
erf.box1.in_box_lo = -3200.   0.     0.
erf.box1.in_box_hi =   800. 400.  1600.

# PROBLEM PARAMETERS
prob.T_0   = 300.0
prob.U_0   = 10.0
prob.rho_0 = 1.16

prob.custom_terrain_type = "WoA"
prob.hmax                = 400.0
prob.L                   = 1000.0
