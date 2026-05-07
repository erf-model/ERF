# ------------------  INPUTS TO MAIN PROGRAM  -------------------
erf.prob_name = "RICO"
stop_time = 3600
max_step = 0

amrex.fpe_trap_invalid = 1
erf.fix_random_seed = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  12800     12800    4000
amr.n_cell           =     32        32      25

geometry.is_periodic = 1 1 0

# MOST BOUNDARY (DEFAULT IS ADIABATIC FOR THETA)

#RICO options
zlo.type = "surface_layer"
erf.surface_layer.flux_type = "rico" # new MOST type
erf.most.ustar  = 0.001229 # Cm
erf.most.tstar  = 0.001094 # Ch
erf.most.qstar  = 0.001133 # Cq
erf.most.z0     = 0.1 # from BOMEX
erf.most.zref   =160.0 # from BOMEX

erf.most.rico.theta_z0 = 300.0
erf.most.rico.qsat_z0 = 0.022

# NOTE: This should have a qv grad too (use hoextrapcc?!)
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt           = 2.5 # fixed time step depending on grid resolution
erf.fixed_mri_dt_ratio = 4

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
erf.plot_int_1      = 2000    # number of timesteps between plotfiles
erf.plot_vars_1     = density \
                      rhotheta \
                      x_velocity \
                      y_velocity \
                      z_velocity \
                      pressure \
                      temp \
                      theta \
                      qt \
                      qp \
                      qv \
                      qc \
                      qsat \
                      qrain \
                      rain_accum \
                      rel_humidity \
                      super_droplets_moisture_number_density \
                      super_droplets_moisture_mass_flux_x \
                      super_droplets_moisture_mass_flux_y \
                      super_droplets_moisture_mass_flux_z \
                      super_droplets_moisture_number_density \
                      super_droplets_moisture_mass_density \
                      super_droplets_moisture_radius \
                      super_droplets_moisture_mass_density_H2O \
                      super_droplets_moisture_mass_flux_x_H2O \
                      super_droplets_moisture_mass_flux_y_H2O \
                      super_droplets_moisture_mass_flux_z_H2O \
                      super_droplets_moisture_aerosol_mass_density_NH4HSO4 \
                      super_droplets_moisture_aerosol_mass_flux_x_NH4HSO4 \
                      super_droplets_moisture_aerosol_mass_flux_y_NH4HSO4 \
                      super_droplets_moisture_aerosol_mass_flux_z_NH4HSO4
particles.disable_plt = true

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = true

erf.use_coriolis    = true
erf.coriolis_3d     = false
erf.latitude        = 14.982176712702886  # f = 0.376e-4 1/s

erf.dycore_horiz_adv_type    = Upwind_3rd
erf.dycore_vert_adv_type     = Upwind_3rd
erf.dryscal_horiz_adv_type   = Upwind_3rd
erf.dryscal_vert_adv_type    = Upwind_3rd
erf.moistscal_horiz_adv_type = Upwind_3rd
erf.moistscal_vert_adv_type  = Upwind_3rd

erf.moisture_model  = "SuperDroplets"
erf.buoyancy_type   = 1

#sdm parameters
super_droplets_moisture.stable_redistribute = true
super_droplets_moisture.place_randomly_in_cells = false

super_droplets_moisture.distribution_type = "uniform"
super_droplets_moisture.diagnostics_interval = 1
super_droplets_moisture.coalescence_kernel = "Halls"
super_droplets_moisture.aerosols = NH4HSO4
super_droplets_moisture.density_scaling = true

super_droplets_moisture.num_initializations = 2

super_droplets_moisture.0.initial_aerosol_distribution_type_NH4HSO4 = "radius_log_normal"
super_droplets_moisture.0.initial_aerosol_mean_radius_NH4HSO4 = 30.0e-9 #kg (0.03 mu-m)
super_droplets_moisture.0.initial_aerosol_std_radius_NH4HSO4 = 0.2468600779315258
super_droplets_moisture.0.initial_number_density = 9.0e7 #m^{-3}
super_droplets_moisture.0.initial_particles_per_cell = 4

super_droplets_moisture.1.initial_aerosol_distribution_type_NH4HSO4 = "radius_log_normal"
super_droplets_moisture.1.initial_aerosol_mean_radius_NH4HSO4 = 140.0e-9 #kg (0.14 mu-m)
super_droplets_moisture.1.initial_aerosol_std_radius_NH4HSO4 = 0.5596157879354227
super_droplets_moisture.1.initial_number_density = 1.5e7 #m^{-3}
super_droplets_moisture.1.initial_particles_per_cell = 4

erf.molec_diff_type = "None"

erf.les_type        = "Smagorinsky"
erf.Cs              = 0.17

#erf.les_type = "Deardorff"
#erf.Ck       = 0.1
#erf.sigma_k  = 1.0
#erf.Ce       = 0.1

erf.Pr_t      = 0.33333333333333
erf.Sc_t      = 0.33333333333333

erf.init_type = "input_sounding"

#erf.restart = chk160000

erf.add_custom_rhotheta_forcing        = true
erf.add_custom_moisture_forcing        = true
erf.add_custom_geostrophic_profile     = true
erf.add_custom_w_subsidence            = true
erf.custom_forcing_uses_primitive_vars = true

# Higher values of perturbations lead to instability
# Instability seems to be coming from BC
prob.U_0_Pert_Mag = 0.01
prob.V_0_Pert_Mag = 0.01
prob.W_0_Pert_Mag = 0.0

prob.pert_ref_height = 1600.0
prob.T_0_Pert_Mag    = 0.1
prob.qv_0_Pert_Mag   = 0.000025

prob.advection_heating_rate   = -2.8935E-5
prob.source_cutoff            = 1500.0
prob.source_cutoff_transition = 1500.0

prob.advection_moisture_rate            = -1.157E-8
prob.moisture_source_cutoff            	= 2980.0
prob.moisture_source_cutoff_transition 	= 200.0

prob.wbar_sub_max    = -0.005
prob.wbar_cutoff_max = 2260.0
prob.wbar_cutoff_min = 2500.0

prob.custom_TKE      = true
