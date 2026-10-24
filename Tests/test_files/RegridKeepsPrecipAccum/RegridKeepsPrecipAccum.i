# Kessler rain from a supersaturated layer below 150 m (sounding_moist), horizontally uniform,
# in a laterally periodic box, so every column rains the same.  A regrid at t = 30 s moves the
# fine grids.  A regrid that keeps the surface rain accumulator leaves it equal on both levels;
# one that rebuilds it starts the newly placed fine columns again from zero.
max_step = 120
erf.prob_name = "ABL"
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_moist"

geometry.prob_extent = 1600. 1600. 1024.
geometry.is_periodic = 1 1 0
amr.n_cell           = 16 16 32
zlo.type = "SlipWall"
zhi.type = "SlipWall"

amr.max_level        = 1
amr.ref_ratio        = 2
erf.dt_ref_ratio     = 1
erf.regrid_int       = 1

erf.refinement_indicators = box1 box2
erf.box1.in_box_lo = 400.  400.    0.
erf.box1.in_box_hi = 800.  800.  512.
erf.box1.end_time  = 30.
erf.box2.in_box_lo = 600.  600.    0.
erf.box2.in_box_hi = 1200. 1200. 512.
erf.box2.start_time = 30.

erf.use_gravity     = true
erf.moisture_model  = "Kessler"
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.pbl_type        = "None"
erf.fixed_dt           = 0.5
erf.fixed_mri_dt_ratio = 6

erf.check_int     = -1
erf.plot_file_1   = plt
erf.plot_int_1    = 120
erf.plot_vars_1   = density theta qv qc qp
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1  = 120
erf.plot2d_vars_1 = precip_rain_accum
erf.v = 0
