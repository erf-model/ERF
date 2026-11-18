# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# A regrid that MOVES an existing fine level must keep the surface that level
# has evolved. ERF::RemakeLevel rebuilds the level on its new grids, and
# init_stuff -> define_level reallocates the two-stream SEB's t_sfc/q_sfc at the
# erf.rad_t_sfc scalar. RemakeLevel therefore retains the old fields, fills the
# new ones from the parent (fill_seb_from_coarse, for cells the old grids never
# covered) and copies the retained values on top.
#
# The other multi-level SEB cases tag a fixed region, so the fine BoxArray never
# changes and RemakeLevel never runs. Here two refinement boxes take turns:
# boxa (coarse x-cells 0-3) until t = 4.75 s, then boxb (coarse x-cells 2-5).
# With erf.regrid_int = 1 the regrid at the start of step 11 moves level 1 from
# fine x-cells 0-7 to 4-11: it keeps 4-7, adds 8-11 and drops 0-3.
#
# check_two_stream_seb_parity.py --regridded-from compares plt2d00011 (new grids)
# with plt2d00009/plt2d00010 (old grids) and plt2d00000:
#   - kept cells: each coarse block's mean, and each fine cell's deviation from
#     it, must extrapolate one step past step 10. Dropping the whole restore
#     resets the mean to erf.rad_t_sfc; dropping only the copy keeps the mean
#     (the interpolation is conservative) but loses the deviations.
#   - added cells: each block's mean must extrapolate the PARENT's surface one
#     step past step 10, as for a newly created level.
#
# Everything else is TwoStream_PrognosticSEBLateLevel's deck (itself
# TwoStream_PrognosticSEBMultiLevel's), on a domain twice as long in x.

erf.prob_name = "ABL"

max_step = 11
stop_time = 10.0
amrex.fpe_trap_invalid = 0

# Deliberately left on AMReX's default MFIter tile size, which splits the
# domain in z. The column sweep must be independent of that tiling; an
# earlier version restarted the sweep at the bottom of every z tile and
# this case is what catches that.

geometry.prob_extent = 2048 1024 1024
amr.n_cell           = 8 4 32
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
# Fine grids that start and stop on coarse-cell boundaries, so each coarse block
# is either wholly kept or wholly added by the regrid (the checker requires it).
amr.blocking_factor = 2

# refine_whole_domain_dir = 2 makes AMReX cluster in the horizontal only and
# emit boxes that span z, which is what the column sweep requires.
amr.refine_whole_domain_dir = 2

# The time windows hand over at t = 4.75 s, half a step before the regrid at
# the start of step 11 (t = 5.0 s), so the handover does not depend on how
# t = 5.0 rounds. Both boxes span z (no in_box z bounds), as the column sweep
# requires.
erf.refinement_indicators = boxa boxb
erf.boxa.max_level = 1
erf.boxa.in_box_lo = 0.0 0.0
erf.boxa.in_box_hi = 1024.0 1024.0
erf.boxa.end_time = 4.75
erf.boxb.max_level = 1
erf.boxb.in_box_lo = 512.0 0.0
erf.boxb.in_box_hi = 1536.0 1024.0
erf.boxb.start_time = 4.75
erf.regrid_int = 1

erf.check_file = chk
erf.check_int = -1

erf.plot_file_1 = plt
erf.plot_int_1 = 11
erf.plot_vars_1 = density theta qsrc_sw qsrc_lw

# The surface state itself, which is what the check reads: every step, so it
# has plt2d00009 and plt2d00010 (old grids) and plt2d00011 (new grids).
erf.plot2d_file_1 = plt2d
erf.plot2d_int_1 = 1
erf.plot2d_vars_1 = seb_t_sfc seb_q_sfc

# Terrain: a Witch of Agnesi hill on a fitted mesh, so every column has its own
# layer thicknesses and the fine level resolves the hill the coarse one cannot.
# Without this the problem is horizontally uniform, both levels compute identical
# surface fluxes and evolve identically -- and the test would pass whether or not
# the levels are actually kept in step, which is worth nothing.
#
# Taller and narrower than the MultiLevel case's hill (hmax 100, L 300), and
# centred at x = 1024 m, on the boundary between the kept and added fine cells:
# the kept cells then carry sub-coarse structure (up to ~4e-4 K from their
# block means) that a restore which only re-interpolated would lose. With the
# gentler hill that loss was barely 3x the check's tolerance; here it is ~6x.
erf.terrain_type         = StaticFittedMesh
erf.terrain_smoothing    = 0
prob.custom_terrain_type = "WoA"
prob.dir                 = 0
prob.hmax                = 200.0
prob.L                   = 150.0

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
