# ------------------  INPUTS TO MAIN PROGRAM  -------------------
erf.prob_name = "MultiSpecies Bubble"

max_step  = 1
stop_time = 3600.0

erf.init_type = MoistBaseState

amrex.fpe_trap_invalid = 1
erf.fix_random_seed = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 20000.0 400.0  10000.0
amr.n_cell           = 100     4      50
geometry.is_periodic = 0 1 0
xlo.type = "SlipWall"
xhi.type = "SlipWall"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 0.5
erf.fixed_mri_dt_ratio = 4

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = -1        # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 100        # number of timesteps between plotfiles
erf.plot_vars_1     = density \
                      rhotheta \
                      rhoQ1 \
                      rhoQ2 \
                      rhoadv_0 \
                      x_velocity \
                      y_velocity \
                      z_velocity \
                      pressure \
                      theta \
                      scalar \
                      temp \
                      pres_hse \
                      dens_hse \
                      pert_pres \
                      pert_dens \
                      eq_pot_temp \
                      qt \
                      qv \
                      qc \
                      qrain \
                      rel_humidity \
                      rain_accum \
                      accum_NaCl \
                      qv_water \
                      qc_water \
                      qt_water \
                      sat_ratio_water \
                      accum_water \
                      qv_agua \
                      qc_agua \
                      qt_agua \
                      sat_ratio_agua \
                      accum_agua \
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
                      super_droplets_moisture_aerosol_mass_density_NaCl \
                      super_droplets_moisture_aerosol_mass_flux_x_NaCl \
                      super_droplets_moisture_aerosol_mass_flux_y_NaCl \
                      super_droplets_moisture_aerosol_mass_flux_z_NaCl \
                      super_droplets_moisture_species_mass_density_H2O \
                      super_droplets_moisture_species_mass_flux_x_H2O \
                      super_droplets_moisture_species_mass_flux_y_H2O \
                      super_droplets_moisture_species_mass_flux_z_H2O \
                      super_droplets_moisture_species_mass_density_water \
                      super_droplets_moisture_species_mass_flux_x_water \
                      super_droplets_moisture_species_mass_flux_y_water \
                      super_droplets_moisture_species_mass_flux_z_water \
                      super_droplets_moisture_species_mass_density_agua \
                      super_droplets_moisture_species_mass_flux_x_agua \
                      super_droplets_moisture_species_mass_flux_y_agua \
                      super_droplets_moisture_species_mass_flux_z_agua
particles.disable_plt = true

# SOLVER CHOICES
erf.use_gravity          = true

erf.dycore_horiz_adv_type    = "Upwind_3rd"
erf.dycore_vert_adv_type     = "Upwind_3rd"
erf.dryscal_horiz_adv_type   = "Upwind_3rd"
erf.dryscal_vert_adv_type    = "Upwind_3rd"
erf.moistscal_horiz_adv_type = "Upwind_3rd"
erf.moistscal_vert_adv_type  = "Upwind_3rd"

# PHYSICS OPTIONS
erf.les_type        = "None"
erf.pbl_type        = "None"
erf.moisture_model  = "SuperDroplets"
erf.buoyancy_type   = 1

erf.molec_diff_type  = "ConstantAlpha"
erf.rho0_trans       = 1.0 # [kg/m^3], used to convert input diffusivities
erf.alpha_T          = 0.0 # [m^2/s]
erf.alpha_C          = 0.0

#sdm parameters
super_droplets_moisture.stable_redistribute = true
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.initial_distribution_type = "uniform"
super_droplets_moisture.diagnostics_interval = 1
super_droplets_moisture.coalescence_kernel = "Halls"
super_droplets_moisture.species = water agua
super_droplets_moisture.aerosols = NaCl
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = "mass_constant"
#super_droplets_moisture.initial_aerosol_distribution_type_NaCl = "mass_exponential"
#super_droplets_moisture.initial_aerosol_min_mass_NaCl = 1.0e-22 #kg
super_droplets_moisture.initial_aerosol_mean_mass_NaCl = 1.0e-19 #kg
super_droplets_moisture.initial_number_density = 1.0e7 #m^{-3}
super_droplets_moisture.initial_particles_per_cell = 1

# PROBLEM PARAMETERS (optional)
# warm bubble input
prob.x_c    = 10000.0
prob.z_c    =  2000.0
prob.x_r    =  2000.0
prob.z_r    =  2000.0
prob.T_0    =   300.0

prob.do_moist_bubble = true
prob.theta_pert  = 2.0
prob.qt_init     = 0.02
prob.eq_pot_temp = 320.0

prob.qv_init_water = 0.0025
prob.qv_init_agua  = 0.01
