# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Rayleigh damping toward an input sounding on a level that a regrid remakes.
#
# A linear sounding (theta = 300 + 0.01 z, U = 5 + 0.01 z, V = 2 + 0.005 z) with Rayleigh
# damping of U, V and T in the top 600 m.  Level 1 starts at 125 < x < 375 m (indicator
# "early", until t = 10) and moves to 375 < x < 750 m (indicator "late", from t = 10), so the
# regrid after t = 10 remakes it through RemakeLevel.  The damping targets of a remade level
# used to be reset to the problem defaults (theta = 300 K, U = V = 0), and the remade level was
# relaxed toward those, well below the sounding's 309 to 310 K in the top 100 m by t = 20
# (erf-model/ERF#5142).  The CTest entry bounds theta in the column through the domain center,
# which lies in the late box, from 900 m up.

erf.prob_name = "ABL"

max_step = 100

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 1000  1000  1000
amr.n_cell           =   16    16    32
amr.max_grid_size    = 32
amr.blocking_factor  = 4
geometry.is_periodic = 1 1 0

# REFINEMENT: one full-height box that moves at t = 10
amr.max_level        = 1
amr.ref_ratio        = 2
erf.regrid_int       = 2
erf.refinement_indicators = early late
erf.early.max_level  = 1
erf.early.in_box_lo  = 125.   250.     0.
erf.early.in_box_hi  = 375.   750.  1000.
erf.early.end_time   = 10.0
erf.late.max_level   = 1
erf.late.in_box_lo   = 375.   250.     0.
erf.late.in_box_hi   = 750.   750.  1000.
erf.late.start_time  = 10.0

# BOUNDARIES
zlo.type = "SlipWall"
zhi.type = "SlipWall"

# INITIALIZATION
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_linear"

# RAYLEIGH DAMPING
erf.rayleigh_damp_U   = true
erf.rayleigh_damp_V   = true
erf.rayleigh_damp_T   = true
erf.rayleigh_dampcoef = 0.2
erf.rayleigh_zdamp    = 600.0

# PHYSICS
erf.use_gravity = true
erf.les_type    = "None"
erf.molec_diff_type = "None"

# TIME STEP CONTROL
erf.fixed_dt           = 0.2
erf.fixed_mri_dt_ratio = 4

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.v            = 1
amr.v            = 1

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1   = plt
erf.plot_int_1    = 100
erf.plot_vars_1   = density x_velocity y_velocity z_velocity theta
