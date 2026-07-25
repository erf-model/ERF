# Production selection regression: request the new field before surf_pres and
# verify the writer restores canonical catalog order in the 2D header.
erf.prob_name = "Taylor-Green Vortex"
erf.init_type = Uniform

max_step = 0

geometry.is_periodic = 1 1 0
geometry.prob_extent = 6.283185307179586 6.283185307179586 6.283185307179586
amr.n_cell = 16 16 16
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.fixed_dt = 0.16
erf.v = 1
amr.v = 1
amr.max_level = 0

erf.check_file = chk
erf.check_int = -1

erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 1
erf.plot2d_vars_1 = sea_level_pressure surf_pres

erf.use_gravity = false
erf.les_type = "None"
erf.molec_diff_type = "None"

prob.rho_0 = 1.0
prob.V_0 = 1.0
