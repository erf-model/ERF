# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# The RANS wall distance on levels that a regrid makes (#4691).
#
# The k-equation RANS deck of PBLH_FinePatch, but its level-1 patch does not exist at the
# start: patch A starts at t = 0.5 s, so the regrid after step 1 adds level 1 from level 0
# (MakeNewLevelFromCoarse); at t = 2.5 s patch A ends and patch B, elsewhere, begins, so the
# regrid after step 3 moves level 1 (RemakeLevel).  Each makes level 1 on new grids, which
# resets walldist to bogus_large_value.  The CTest entry bounds walldist in the final
# plotfile over both levels: half a cell (8 m on level 1) to the top cell centre of level 0;
# RANS_WallDist_NewLevel.i stops after step 2 and checks the level added from level 0.

erf.prob_name = "ABL"

max_step = 4

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 3200  3200  1024
amr.n_cell           =   32    32    32
amr.max_grid_size    = 32
amr.blocking_factor  = 4
geometry.is_periodic = 1 1 0

# REFINEMENT: level 1 is added after step 1 and moved after step 3
amr.max_level        = 1
amr.ref_ratio        = 2
erf.regrid_int       = 1
erf.refinement_indicators = patchA patchB
erf.patchA.max_level  = 1
erf.patchA.in_box_lo  =  800.  800.   0.
erf.patchA.in_box_hi  = 1600. 1600. 256.
erf.patchA.start_time = 0.5
erf.patchA.end_time   = 2.5
erf.patchB.max_level  = 1
erf.patchB.in_box_lo  = 1600. 1600.   0.
erf.patchB.in_box_hi  = 2400. 2400. 256.
erf.patchB.start_time = 2.5

# SURFACE LAYER
zlo.type = "surface_layer"
erf.most.z0             = 0.1
erf.most.average_policy = 1

zhi.type       = "SlipWall"
zhi.theta_grad = 0.003

# INITIALIZATION
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_inversion"
prob.KE_0            = 0.5
prob.KE_decay_height = 500.0

# PHYSICS
erf.use_gravity = true

# TURBULENCE
erf.les_type  = "None"
erf.rans_type = "kEqn"
erf.rans_lscale_min = 1.0
erf.max_geom_lscale = 1000.0

# TIME STEP CONTROL
erf.fixed_dt           = 1.0
erf.fixed_mri_dt_ratio = 6

# DIAGNOSTICS & VERBOSITY
erf.v = 1
amr.v = 1

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 4
erf.plot_vars_1 = density theta KE walldist
