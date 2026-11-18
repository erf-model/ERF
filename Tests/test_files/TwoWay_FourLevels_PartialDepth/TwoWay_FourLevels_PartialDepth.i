# Four levels (amr.ref_ratio_vect = 2 2 1 on each), each refined level stopping below the domain
# top, anelastic with two-way coupling, in a uniform flow (erf-model/ERF#4988).
#
# With amr.iterate_grids (the default) AmrMesh::MakeNewGrids builds the levels one at a time and
# then regenerates them all from level 0 once the finest level exists.  Levels 1 and 2 grow by the
# proper-nesting buffer in that pass and are made from scratch a second time; level 3 sits well
# inside level 2 and is not.  The run log shows it (erf.v = 1):
#   BA FROM SCRATCH AT LEVEL 2 ...   (twice)
#   BA FROM SCRATCH AT LEVEL 3 ...   (once)
# When level 2 is made again it must rebuild the flux register on its interface with level 3,
# which was built from its first grids, and fill its own ghost cells, which the second tagging pass
# otherwise fills from a level 1 whose ghost cells are still empty.  The flow is uniform, so every
# step must leave it as it started on every level (plt00000 against plt00003).
max_step  = 3
erf.fixed_dt = 1.6

geometry.prob_lo     = -4800. -3240.    0.
geometry.prob_hi     =  4800.  3360. 1600.
amr.n_cell           =    40     28   40
amr.max_level        = 3
amr.ref_ratio_vect   = 2 2 1  2 2 1  2 2 1
amr.n_error_buf      = 0
amr.blocking_factor  = 2
amr.max_grid_size    = 16

erf.refinement_indicators = box1 box2 box3
erf.box1.max_level = 1
erf.box1.in_box_lo_indices_crse =  11   5  0
erf.box1.in_box_hi_indices_crse =  25  21 27
erf.box2.max_level = 2
erf.box2.in_box_lo_indices_crse =  37  17  0
erf.box2.in_box_hi_indices_crse =  62  36 25
erf.box3.max_level = 3
erf.box3.in_box_lo_indices_crse =  82  38  0
erf.box3.in_box_hi_indices_crse = 117  69 17

erf.coupling_type = TwoWay

geometry.is_periodic = 1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.compressibility_type = Anelastic
erf.vert_implicit  = true
erf.timestepping_type = MidPoint
erf.use_gravity    = true

erf.molec_diff_type = "None"
erf.les_type  = "None"

erf.init_type = "uniform"
erf.prob_name = "ABL"
prob.rho_0 = 1.225
prob.T_0   = 300.0
prob.U_0   = 10.0
erf.theta_ref = 300.0

erf.sum_interval = -1
erf.v            = 1   # prints BA FROM SCRATCH AT LEVEL <lev>, which the test counts
erf.plot_file_1  = plt
erf.plot_int_1   = 3
erf.plot_vars_1  = density rhotheta x_velocity y_velocity z_velocity
erf.check_int    = -1
