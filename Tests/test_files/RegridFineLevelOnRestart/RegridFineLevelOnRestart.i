# Restart parity across a regrid that re-makes level 1 over the same cells on a new box
# layout (erf-model/ERF#5228, #5236).  Level 1 covers x < 500 until t = 10 and x > 500 after;
# the restart leg changes amr.max_grid_size (or the rank count), so the first regrid after
# the restart re-chops level 1.  The re-made level must get the same base state a level built
# from scratch on those grids has, or the answer changes.
erf.prob_name = "Regrid of a fine level on restart"

max_step = 40

geometry.prob_lo = 0 0 0
geometry.prob_hi = 1000 1000 1000
amr.n_cell       = 16 16 32
geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"

erf.fixed_dt           = 0.5
erf.fixed_mri_dt_ratio = 8

amr.max_level       = 1
amr.ref_ratio       = 2
erf.regrid_int      = 2
amr.max_grid_size   = 16
amr.blocking_factor = 8

erf.refinement_indicators = boxA boxB
erf.boxA.in_box_lo  =   0 0    0
erf.boxA.in_box_hi  = 500 1000 1000
erf.boxA.end_time   = 10.0
erf.boxB.in_box_lo  = 500 0    0
erf.boxB.in_box_hi  = 1000 1000 1000
erf.boxB.start_time = 10.0

erf.check_file  = chk
erf.plot_file_1 = plt
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta

erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.use_gravity     = true

erf.v = 1
amr.v = 1
