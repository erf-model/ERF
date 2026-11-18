# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Custom rhotheta, moisture and subsidence forcing on a refined level that a regrid creates and
# later remakes on larger grids.
#
# A cut-down Exec/CanonicalTests/Bomex deck.  Level 1 does not exist at start-up: the "first"
# box appears at t = 0.75, so the regrid after step 3 (t = 0.9) makes level 1 with
# MakeNewLevelFromCoarse; the "second" box appears at t = 1.95, so the regrid after step 7
# (t = 2.1) remakes it with RemakeLevel on twice as many cells and boxes.  The forcing was built only in
# InitData_post, so the new level had null sources and an empty subsidence profile, and the
# remade level kept its sources on the old grids (erf-model/ERF#5151).  advance_dycore now
# aborts if the forcing is not built on the level's current grids; the CTest entry runs the
# deck and bounds theta.

erf.prob_name = "BOMEX"

max_step = 10

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 6400 6400 4000
amr.n_cell           =   16   16   32
amr.max_grid_size    = 16 8   # level 1 goes from 1 to 2 boxes when it is remade
amr.blocking_factor  = 4
geometry.is_periodic = 1 1 0

# REFINEMENT: level 1 made at t = 0.9 and grown at t = 2.1
amr.max_level        = 1
amr.ref_ratio        = 2
erf.regrid_int       = 1
erf.refinement_indicators = first second
erf.first.max_level   = 1
erf.first.in_box_lo   = 1600. 1600.    0.
erf.first.in_box_hi   = 3200. 3200. 1000.
erf.first.start_time  = 0.75
erf.second.max_level  = 1
erf.second.in_box_lo  = 1600. 1600.    0.
erf.second.in_box_hi  = 4800. 3200. 1000.
erf.second.start_time = 1.95

# MOST BOUNDARY
zlo.type = "surface_layer"
erf.surface_layer.flux_type = "custom"
erf.most.ustar  = 0.28
erf.most.tstar  = 8.0e-3
erf.most.qstar  = 5.2e-5
erf.most.z0     = 0.1
erf.most.zref   = 100.0   # at least dz/2 = 62.5 on this grid

zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt           = 0.3
erf.fixed_mri_dt_ratio = 4

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.v            = 1
amr.v            = 1

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 10
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta qv

erf.use_gravity  = true
erf.use_coriolis = true
erf.coriolis_3d  = false
erf.latitude     = 14.982176712702886

erf.dycore_horiz_adv_type    = Upwind_3rd
erf.dycore_vert_adv_type     = Upwind_3rd
erf.dryscal_horiz_adv_type   = Upwind_3rd_SL
erf.dryscal_vert_adv_type    = Upwind_3rd_SL
erf.moistscal_horiz_adv_type = Upwind_3rd_SL
erf.moistscal_vert_adv_type  = Upwind_3rd_SL

erf.moisture_model = "SAM"
erf.buoyancy_type  = 1

erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.17
erf.Pr_t            = 0.33333333333333
erf.Sc_t            = 0.33333333333333

erf.init_type     = "input_sounding"
erf.sounding_type = Ideal

erf.add_custom_rhotheta_forcing        = true
erf.add_custom_moisture_forcing        = true
erf.add_custom_geostrophic_profile     = true
erf.add_custom_w_subsidence            = true
erf.custom_forcing_uses_primitive_vars = true

prob.U_0_Pert_Mag = 0.01
prob.V_0_Pert_Mag = 0.01
prob.W_0_Pert_Mag = 0.0

prob.pert_ref_height = 1600.0
prob.T_0_Pert_Mag    = 0.1
prob.qv_0_Pert_Mag   = 0.000025

prob.advection_heating_rate   = -2.3148E-5
prob.source_cutoff            = 1500.0
prob.source_cutoff_transition = 1500.0

prob.advection_moisture_rate           = -1.2E-8
prob.moisture_source_cutoff            = 300.0
prob.moisture_source_cutoff_transition = 200.0

prob.wbar_sub_max    = -0.0065
prob.wbar_cutoff_max = 1500.0
prob.wbar_cutoff_min = 2100.0

prob.custom_TKE = true
