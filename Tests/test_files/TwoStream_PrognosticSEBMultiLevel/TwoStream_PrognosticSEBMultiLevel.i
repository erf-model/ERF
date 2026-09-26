# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# The prognostic surface energy balance on a refined hierarchy. Every level
# evolves its own force-restore surface temperature, and that temperature is the
# longwave boundary condition -- so the levels must agree about the ground they
# share. ERF averages t_sfc and q_sfc down after the finer levels advance.
#
# The refinement criterion tags EVERY cell (value_less is above every theta in
# the column), so level 1 covers the whole domain and every coarse cell is
# covered by exactly ref_ratio^2 fine cells.
#
# check_two_stream_seb_parity.py asserts the invariant that average_down
# establishes exactly: each coarse seb_t_sfc equals the mean of the fine cells
# above it. That is an equality to round-off when the average-down runs, so the
# check uses a round-off tolerance rather than a physical one. Without the
# average-down level 0 keeps its own coarse-grid solution, which differs from
# the fine mean by the discretization error -- orders of magnitude above the
# tolerance, and detectable from the first output. Comparing domain MEANS would
# not work here: average_down is mean-preserving, so the means agree either way.

erf.prob_name = "ABL"

max_step = 6
stop_time = 10.0
amrex.fpe_trap_invalid = 0

# Deliberately left on AMReX's default MFIter tile size, which splits the
# domain in z. The column sweep must be independent of that tiling; an
# earlier version restarted the sweep at the bottom of every z tile and
# this case is what catches that.

geometry.prob_extent = 1024 1024 1024
amr.n_cell           = 4 4 32
amr.max_grid_size_z = 128
geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"
zhi.theta_grad = 0.003

erf.fixed_dt = 0.5
erf.sum_interval = 1
erf.v = 1
amr.v = 1
amr.max_level = 1
amr.ref_ratio_vect = 2 2 1
amr.n_error_buf = 0

# Tag on theta, which increases with height (zhi.theta_grad above): the tagged
# region is the lower part of the column and would give a patch that stops
# short of the domain top. refine_whole_domain_dir = 2 makes AMReX cluster in
# the horizontal only and emit boxes that span z, which is what the column
# sweep requires (and what the model's abort message recommends).
amr.refine_whole_domain_dir = 2
erf.refinement_indicators = lowth
erf.lowth.max_level = 1
erf.lowth.field_name = theta
erf.lowth.value_less = 1.0e4   # above every theta: tags the whole domain

erf.check_file = chk
erf.check_int = -1

erf.plot_file_1 = plt
erf.plot_int_1 = 2
erf.plot_vars_1 = density theta qsrc_sw qsrc_lw

# The surface state itself, which is what the parity check reads.
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 2
erf.plot2d_vars_1 = seb_t_sfc seb_q_sfc

# Terrain: a Witch of Agnesi hill on a fitted mesh, so every column has its own
# layer thicknesses and the fine level resolves the hill the coarse one cannot.
# Without this the problem is horizontally uniform, both levels compute identical
# surface fluxes and evolve identically -- and the test would pass whether or not
# the levels are actually kept in step, which is worth nothing.
erf.terrain_type         = StaticFittedMesh
erf.terrain_smoothing    = 0
prob.custom_terrain_type = "WoA"
prob.dir                 = 0
prob.hmax                = 100.0
prob.L                   = 300.0

erf.use_gravity = true
erf.molec_diff_type = "None"
erf.les_type = "None"
erf.pbl_type = "None"
erf.theta_ref = 300.0

erf.init_type = "input_sounding"
erf.sounding_type = Ideal
erf.input_sounding_file = "input_sounding"

erf.use_coriolis = false
erf.abl_driver_type = "None"

# RADIATION - TwoStream, SW + LW, clear sky, fixed sun
erf.radiation_model = "TwoStream"
erf.radiation.sw_enabled = true
erf.radiation.lw_enabled = true
erf.radiation.tau_per_layer = 0.00625
erf.radiation.tau_lw_per_layer = 1.0
erf.fixed_solar_zenith_angle = 0.5    # cos(60 deg): the cosine, as RRTMGP takes it
erf.fixed_total_solar_irradiance = 1361.0
erf.rad_t_sfc = 300.0    # surface temperature [K] where no LSM or surface layer supplies one (shared with RRTMGP)
erf.radiation.v = 0
erf.radiation.diag_csv_enable = true
erf.radiation.diag_stdout_enable = true

# The prognostic surface energy balance: the feature under test.
erf.radiation.seb_enable = true
erf.radiation.seb_prognostic_enable = true
erf.radiation.diag_enable = true

# The surface has to actually move, or the test passes on a constant and proves
# nothing: with the defaults every flux in the balance is zero and T_s sits at
# erf.rad_t_sfc forever. Drive it from the sweep's own surface fluxes and shorten
# the restore timescale so the drift is resolvable in a few steps.
erf.radiation.seb_use_radiation_fluxes = true
erf.radiation.seb_restore_timescale_s = 3600.0
