# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# The TwoStream_ColumnHeating column over a Witch of Agnesi hill. The layer
# thicknesses come from the interface heights of the terrain-fitted mesh, so
# this deck covers the non-uniform-grid path, and because the columns differ
# the runner's comparison of the diagnostics CSV between 1 rank and NRANKS
# fails if any domain mean is formed without an MPI reduction.
#
# The TwoStream radiation sweep must treat k = 0 as the surface layer and
# the highest k as the top of the atmosphere, and the LW heating sign must
# give cooling to space. This short SW + LW run writes the per-level heating
# rates (qsrc_sw, qsrc_lw) to plt00002; TwoStreamRadiationCheck verifies the
# vertical structure: SW heating strongest at the top layer, LW cooling
# strongest at the top layer, net LW cooling of the column.
erf.prob_name = "ABL"

max_step = 2
stop_time = 10.0
amrex.fpe_trap_invalid = 0

# Deliberately left on AMReX's default MFIter tile size, which splits the
# domain in z. The column sweep must be independent of that tiling; an
# earlier version restarted the sweep at the bottom of every z tile and
# this case is what catches that.

geometry.prob_lo     = -768 -512    0
geometry.prob_hi     =  768  512 1024
amr.n_cell           = 6 4 32
# Two boxes of unequal width (4 and 2 columns) that are not mirror images
# of each other about the hill, so each rank's share of the columns has a
# different mean from the whole domain.
amr.max_grid_size_x = 4
amr.blocking_factor = 2
amr.max_grid_size_z = 128
geometry.is_periodic = 1 1 0

# TERRAIN: axisymmetric Witch of Agnesi hill centred in the domain on a
# terrain-fitted mesh, so every column has its own layer thicknesses (taken
# from the nodal heights) and the domain-mean diagnostics differ from any
# single rank's share of the columns.
erf.terrain_type         = StaticFittedMesh
erf.terrain_smoothing    = 0
prob.custom_terrain_type = "WoA"
prob.dir                 = 2
prob.hmax                = 100.0
prob.L                   = 300.0

zlo.type = "SlipWall"
zhi.type = "SlipWall"
zhi.theta_grad = 0.003

erf.fixed_dt = 0.5
erf.sum_interval = 1
erf.v = 1
amr.v = 1
amr.max_level = 0

erf.check_file = chk
erf.check_int = -1

erf.plot_file_1 = plt
erf.plot_int_1 = 2
erf.plot_vars_1 = density theta qsrc_sw qsrc_lw

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
erf.radiation.diag_csv_enable = false
erf.radiation.diag_stdout_enable = false
