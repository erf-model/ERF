# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# PBL height on a refined level that does not reach the ground.
#
# The deck of PBLH_FinePatch (6 K inversion at 500-600 m, k-equation RANS with the length scale
# capped from the MYNN25 PBL height) with the patch aloft, from 256 m to 768 m: the level has no
# column that starts at the ground to scan, so it must take the height of level 0 (~461 m)
# rather than zero, which would collapse the length cap to erf.rans_lscale_min on the patch.
# The CTest entry bounds the pblh of the 2D plotfile over both levels.

erf.prob_name = "ABL"

max_step = 4

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 3200  3200  1024
amr.n_cell           =   32    32    32
amr.max_grid_size    = 32
amr.blocking_factor  = 4
geometry.is_periodic = 1 1 0

# REFINEMENT: one patch aloft, across the inversion
amr.max_level        = 1
amr.ref_ratio        = 2
erf.refinement_indicators = patch
erf.patch.max_level  = 1
erf.patch.in_box_lo  =  800.  800. 256.
erf.patch.in_box_hi  = 2400. 2400. 768.

# SURFACE LAYER
zlo.type = "surface_layer"
erf.most.z0             = 0.1
erf.most.pblh_calc      = "MYNN25"
erf.most.average_policy = 1

zhi.type       = "SlipWall"
zhi.theta_grad = 0.003

# INITIALIZATION
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_inversion"
# uniform: the KE profile option (prob.KE_decay_height) needs every box to start at the ground
prob.KE_0            = 0.5

# PHYSICS
erf.use_gravity = true

# TURBULENCE
erf.les_type  = "None"
erf.rans_type = "kEqn"
erf.rans_lscale_from_pblh = true
erf.rans_lscale_min       = 1.0
erf.max_geom_lscale       = 1000.0

# TIME STEP CONTROL
erf.fixed_dt           = 1.0
erf.fixed_mri_dt_ratio = 6

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.v            = 1
amr.v            = 1

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1   = plt
erf.plot_int_1    = 4
erf.plot_vars_1   = density x_velocity y_velocity z_velocity theta KE Kmv
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1  = 4
erf.plot2d_vars_1 = pblh u_star
