# ParentProject_HelperRegression: the case of Tests/StructureModelExample/inputs_structure_example
# (one point pushing on still air with 50 N along x), kept here as a parent repository keeps
# its own decks under test_files/<test name>. Only the output differs: one plotfile at the
# last step and no checkpoint. Keep the other lines the same as the model's deck: the gold
# is the plotfile ParentProject_StructureExample writes from that deck.

max_step  = 16
stop_time = 1.0e9

geometry.prob_lo     = 0.0 0.0 0.0
geometry.prob_hi     = 400.0 400.0 400.0
geometry.is_periodic = 1 1 0
amr.n_cell           = 16 16 16
amr.max_level        = 0
amr.max_grid_size    = 8 8 16
amr.blocking_factor  = 2

zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.init_type = "Uniform"   # constant density and temperature, at rest
prob.rho_0    = 1.0
prob.T_0      = 300.0

erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.use_gravity     = false
erf.use_coriolis    = false

erf.fixed_dt     = 0.25   # acoustic substeps stay below c dtau / dx = 0.3
erf.v            = 0
amr.v            = 0
erf.sum_interval = -1

erf.plot_file_1 = plt
erf.plot_int_1  = 16
erf.plot_vars_1 = density x_velocity y_velocity z_velocity structure_example_fx structure_example_fy structure_example_fz
erf.check_file  = chk
erf.check_int   = -1

erf.structures.model = example

structure.example.position = 200.0 200.0 200.0
structure.example.force    = 50.0 0.0 0.0
structure.example.epsilon  = 50.0
