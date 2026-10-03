# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Input sponge on a refined level that does not reach the top of the domain.
#
# A neutral layer at 8 m/s under a 6 K inversion (sounding_inversion), damped toward 6 m/s
# (sponge_profile) over x > 2400 m.  Level 1 is a patch at the ground, inside the sponge region,
# that ends at 256 m.  The sponge reference profile of each level is set from the largest cell
# height at every k of that level; the reduction used to read every k from every box, so the
# patch, whose boxes stop at k = 15 of 64, was read outside its data at start-up (a crash).
# Only the xhi side is switched on and amrex.fpe_trap_invalid = 1: the sponge kernels divided by
# the ramp lengths of the sides that are switched off too (evaluated speculatively by the
# optimiser), whose ends were never set, and the first step trapped even on one level.
# The CTest entry runs it and bounds x_velocity (6 m/s target, 8 m/s free stream plus the start-up
# adjustment, 8.24 m/s).

erf.prob_name = "ABL"

max_step = 4

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 3200  3200  1024
amr.n_cell           =   32    32    32
amr.max_grid_size    = 32
amr.blocking_factor  = 4
geometry.is_periodic = 1 1 0

# REFINEMENT: one patch at the ground inside the sponge region
amr.max_level        = 1
amr.ref_ratio        = 2
erf.refinement_indicators = patch
erf.patch.max_level  = 1
erf.patch.in_box_lo  = 1600.  800.   0.
erf.patch.in_box_hi  = 3200. 2400. 256.

# BOUNDARIES
zlo.type = "SlipWall"
zhi.type = "SlipWall"

# INITIALIZATION
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_inversion"

# SPONGE
erf.sponge_type            = "input_sponge"
erf.input_sponge_file      = "sponge_profile"
erf.sponge_strength        = 0.1
erf.use_xhi_sponge_damping = true
erf.xhi_sponge_start       = 2400.0

# PHYSICS
erf.use_gravity = true
erf.les_type    = "None"
erf.molec_diff_type = "None"

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
erf.plot_vars_1   = density x_velocity y_velocity z_velocity theta
