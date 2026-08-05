# ------------------  INPUTS TO MAIN PROGRAM  -------------------
erf.prob_name = "Taylor-Green Vortex"
erf.init_type = Uniform

max_step = 0

amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.is_periodic = 1 1 0
geometry.prob_extent = 6.283185307179586 6.283185307179586 6.283185307179586
amr.n_cell = 16 16 16
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.fixed_dt = 0.16
erf.fixed_mri_dt_ratio = 6
erf.sum_interval = 1
erf.v = 1
amr.v = 1
amr.max_level = 0

erf.check_file = chk
erf.check_int = -1

erf.plot_file_1 = plt
erf.plot_int_1 = 1
erf.plot_vars_1 = density qsrc_sw walldist u_t_avg Kmv diss theta VPD

erf.use_gravity = false
erf.les_type = "None"
erf.molec_diff_type = "None"

prob.rho_0 = 1.0
prob.V_0 = 1.0
