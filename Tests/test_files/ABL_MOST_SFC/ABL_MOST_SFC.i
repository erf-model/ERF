# ------------------  INPUTS TO MAIN PROGRAM  -------------------
erf.prob_name = "ABL"

max_step = 10

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  3200   3200   4000
amr.n_cell           =    32     32    100

geometry.is_periodic = 1 1 0

# MOST BOUNDARY (DEFAULT IS ADIABATIC FOR THETA)
zlo.type = "surface_layer"

# Surface time-varying flux forcing file
erf.most.use_sfc_fluxes = true
erf.most.sfc_file = "sfc.txt"

erf.surface_layer.flux_type = "custom"
erf.most.ustar  = 0.0
erf.most.tstar  = 0.0
erf.most.qstar  = 0.0
erf.most.z0     = 0.1 # from BOMEX
erf.most.zref   = 20.0 # from BOMEX

erf.is_land = 0

# NOTE: This should have a qv grad too (use hoextrapcc?!)
zhi.type = "SlipWall"

# RAYLEIGH DAMPING
erf.rayleigh_damp_W = true
erf.rayleigh_damp_U = true
erf.rayleigh_damp_V = true
erf.rayleigh_damp_T = true
erf.rayleigh_dampcoef = 0.2
erf.rayleigh_zdamp = 500.

# TIME STEP CONTROL
erf.fixed_dt = 1.0
erf.fixed_fast_dt = 0.2

# large scale forcing
#erf.large_scale_forcing = true
#erf.forcing_timescale = 3600.0 #tauls
#erf.large_scale_forcing_file = "lsf.txt"

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp
erf.data_log       = "surf" "mean" "flux" "subgrid"
erf.profile_int    = 200     # (every minute with dt = 0.075)

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk     # root name of checkpoint file
erf.check_int       = -1      # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt     # prefix of plotfile name
erf.plot_int_1      = 10      # number of timesteps between plotfiles
erf.plot_vars_1     = density rhotheta x_velocity y_velocity z_velocity pressure temp theta qv

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = true

erf.use_coriolis    = false
erf.coriolis_3d     = false

erf.dycore_horiz_adv_type    = Upwind_3rd
erf.dycore_vert_adv_type     = Upwind_3rd
erf.dryscal_horiz_adv_type   = Upwind_3rd_SL
erf.dryscal_vert_adv_type    = Upwind_3rd_SL
erf.moistscal_horiz_adv_type = Upwind_3rd_SL
erf.moistscal_vert_adv_type  = Upwind_3rd_SL

erf.moisture_model  = "Morrison"
erf.buoyancy_type   = 2

erf.molec_diff_type = "None"

erf.les_type        = "Smagorinsky"
erf.Cs              = 0.17

erf.Pr_t      = 0.33333333333333
erf.Sc_t      = 0.33333333333333

erf.init_type = "input_sounding"
erf.input_sounding_file = "input_sounding"

erf.add_custom_rhotheta_forcing        = false
erf.add_custom_moisture_forcing        = false
erf.add_custom_geostrophic_profile     = false
erf.add_custom_w_subsidence            = false
erf.custom_forcing_uses_primitive_vars = false

# Higher values of perturbations lead to instability
# Instability seems to be coming from BC
prob.U_0_Pert_Mag = 0.01
prob.V_0_Pert_Mag = 0.01
prob.W_0_Pert_Mag = 0.0

prob.pert_ref_height = 1600.0
prob.T_0_Pert_Mag    = 0.1
