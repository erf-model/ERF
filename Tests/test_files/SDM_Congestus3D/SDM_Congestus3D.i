# ------------------  INPUTS TO MAIN PROGRAM  -------------------
erf.prob_name = "SDM_Congestus3D"

stop_time = 8500
max_step = 20

amrex.fpe_trap_invalid = 1
erf.fix_random_seed = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent    =   6400.  6400.  10000.
amr.n_cell              =     32     32      50

amrex.async_out = 1

geometry.is_periodic = 1 1 0
zhi.type = "SlipWall"

# MOST BOUNDARY (DEFAULT IS ADIABATIC FOR THETA)
zlo.type = "surface_layer"
erf.surface_layer.flux_type = "custom"
erf.most.ustar  = 0.28   # ustar
erf.most.tstar  = 0.0 # theta flux
erf.most.qstar  = 0.0 # qv    flux
erf.most.z0     = 0.1
erf.most.zref   = 200.0

# INITIALIZATION
erf.init_type = input_sounding

# TIME STEP CONTROL
erf.fixed_dt = 0.5
erf.fixed_fast_dt = 0.1

# DIAGNOSTICS & VERBOSITY
amr.v               = 1     # verbosity in Amr.cpp
erf.v               = 1     # verbosity in ERF.cpp -- needs to be 1 to write out data_log files
erf.sum_interval    = 1     # timesteps between computing mass

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk  # root name of checkpoint file
erf.check_int       = -1     # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt # prefix of plotfile name
erf.plot_int_1      = 120 # number of timesteps between plotfiles (DEBUG)
erf.plot_vars_1     = density \
                      rhotheta \
                      x_velocity \
                      y_velocity \
                      z_velocity \
                      pressure \
                      theta \
                      temp \
                      rhoQ1 \
                      qt \
                      qv \
                      qc \
                      qp \
                      qsat \
                      qrain \
                      rain_accum \
                      accum_NH42SO4 \
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
                      super_droplets_moisture_aerosol_mass_density_NH42SO4 \
                      super_droplets_moisture_aerosol_mass_flux_x_NH42SO4 \
                      super_droplets_moisture_aerosol_mass_flux_y_NH42SO4 \
                      super_droplets_moisture_aerosol_mass_flux_z_NH42SO4
particles.disable_plt = true

# SOLVER CHOICES
erf.use_gravity = true

erf.moisture_model  = "SuperDroplets"
erf.buoyancy_type   = 1

#sdm parameters
super_droplets_moisture.stable_redistribute = true
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.initial_distribution_type = "uniform"
super_droplets_moisture.diagnostics_interval = 1
super_droplets_moisture.coalescence_kernel = "Halls"
super_droplets_moisture.aerosols = NH42SO4
super_droplets_moisture.density_scaling = true

super_droplets_moisture.num_initializations = 1
super_droplets_moisture.initial_aerosol_distribution_type_NH42SO4 = "mass_constant"
super_droplets_moisture.initial_aerosol_mean_mass_NH42SO4 = 2.649201195e-17

erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.17

erf.add_custom_rhotheta_forcing = true
erf.add_custom_moisture_forcing = true
erf.spatial_rhotheta_forcing = true
erf.spatial_moisture_forcing = true
erf.custom_forcing_uses_primitive_vars = true
prob.advection_heating_rate   = 3e-3 #trapp has max of 3e-3 K.kg/(m3s)
prob.advection_moisture_rate  = 1e-6 #trapp has 1e-6 kg/(m3s)

prob.advection_heating_rate_base   = 1e-3 #trapp has max of 3e-3 K.kg/(m3s)
prob.advection_moisture_rate_base  = 4e-7 #trapp has 1e-6 kg/(m3s)

prob.pert_deltaU = 0.1
prob.pert_deltaV = 0.1
prob.pert_deltaT = 0.01
prob.pert_deltaQV = 2.5e-5
prob.pert_periods_T = 10.0
prob.pert_periods_QV = 10.0
