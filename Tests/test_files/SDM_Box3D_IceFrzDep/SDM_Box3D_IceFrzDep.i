# ------------------  INPUTS TO MAIN PROGRAM  -------------------
erf.prob_name = "Bubble"
max_step = 10
stop_time = 5.0

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =  0.   0.   0.
geometry.prob_hi     =  0.04 0.04 0.04
amr.n_cell           =  4    4    4      # dx=dy=dz=0.01m
geometry.is_periodic =  1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt       = 0.00001     # fixed time step [s]
erf.substepping_type = "None"

# initialization type
erf.init_type = "Uniform"

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = -1         # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 1          # number of timesteps between plotfiles
erf.plot_vars_1     = density \
                      rhotheta \
                      rhoQ1 \
                      rhoQ2 \
                      rhoQ3 \
                      rhoQ4 \
                      rhoQ5 \
                      rhoQ6 \
                      qt \
                      qv \
                      qc \
                      qrain \
                      qi \
                      qgraup \
                      qsnow \
                      qsat \
                      rain_accum \
                      x_velocity \
                      y_velocity \
                      z_velocity \
                      pressure \
                      theta \
                      temp \
                      pres_hse \
                      dens_hse \
                      pert_pres \
                      pert_dens \
                      rel_humidity \
                      super_droplets_moisture_number_density \
                      super_droplets_moisture_mass_flux_y \
                      super_droplets_moisture_mass_flux_z \
                      super_droplets_moisture_number_density \
                      super_droplets_moisture_mass_density \
                      super_droplets_moisture_radius \
                      super_droplets_moisture_mass_density_H2O \
                      super_droplets_moisture_mass_flux_x_H2O \
                      super_droplets_moisture_mass_flux_y_H2O \
                      super_droplets_moisture_mass_flux_z_H2O \
                      super_droplets_moisture_species_mass_density_ice \
                      super_droplets_moisture_species_mass_flux_x_ice \
                      super_droplets_moisture_species_mass_flux_y_ice \
                      super_droplets_moisture_species_mass_flux_z_ice \
                      super_droplets_moisture_aerosol_mass_density_NaCl \
                      super_droplets_moisture_aerosol_mass_flux_x_NaCl \
                      super_droplets_moisture_aerosol_mass_flux_y_NaCl \
                      super_droplets_moisture_aerosol_mass_flux_z_NaCl 


# SOLVER CHOICES
erf.use_gravity = false
erf.use_coriolis = false

# PHYSICS OPTIONS
erf.les_type        = "None"
erf.pbl_type        = "None"
erf.moisture_model  = "SuperDroplets"
erf.buoyancy_type   = 1

# Super Droplets Options
super_droplets_moisture.stable_redistribute = true
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.initial_distribution_type = "uniform"
super_droplets_moisture.diagnostics_interval = 1
super_droplets_moisture.advect_with_flow = false
super_droplets_moisture.advect_with_gravity = false
super_droplets_moisture.include_coalescence = false

super_droplets_moisture.aerosols = NaCl
super_droplets_moisture.multiplicity_type = "constant"
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.initial_aerosol_mean_mass_NaCl = 1.0e-19 #kg
super_droplets_moisture.initial_species_distribution_type_H2O = "mass_constant"
super_droplets_moisture.initial_species_mean_mass_H2O = 1.0e-14 #kg
super_droplets_moisture.initial_number_density = 1.0e7 #m^{-3}
super_droplets_moisture.initial_particles_per_cell = 1

# PROBLEM PARAMETERS (optional)
prob.U_0    = 0.0
prob.x_c    = 2.0
prob.z_c    = 2.0
prob.x_r    = 1e99
prob.z_r    = 1e99

prob.do_moist_bubble = true
prob.theta_pert  = -60.0
prob.qt_init     = 0.02
prob.eq_pot_temp = 300.0
