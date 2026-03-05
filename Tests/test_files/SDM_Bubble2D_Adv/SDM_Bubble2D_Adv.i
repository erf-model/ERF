# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step  = 50
stop_time = 3600.0

erf.init_type = Isentropic

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
erf.plot_int_1      = 50         # number of timesteps between plotfiles
erf.plot_vars_1     = density \
                      rhotheta \
                      x_velocity \
                      y_velocity \
                      z_velocity \
                      pressure \
                      theta \
                      scalar \
                      temp \
                      eq_pot_temp \
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
                      super_droplets_moisture_aerosol_mass_flux_z_NaCl

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
erf.dynamicViscosity = 0.0 # [kg/(m-s)] ==> nu = 75.0 m^2/s
erf.alpha_T          = 0.0 # [m^2/s]
erf.alpha_C          = 0.0

#sdm parameters
super_droplets_moisture.stable_redistribute = true
super_droplets_moisture.initial_distribution_type = "uniform"
super_droplets_moisture.diagnostics_interval = 1
super_droplets_moisture.include_phase_change = false
super_droplets_moisture.include_coalescence = false
super_droplets_moisture.advect_with_gravity = false
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.aerosols = NaCl
super_droplets_moisture.multiplicity_type = "constant"
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.initial_aerosol_mean_mass_NaCl = 1.0e-19 #kg
super_droplets_moisture.initial_particles_per_cell = 1

# PROBLEM PARAMETERS (optional)
# warm bubble input
prob.T_pert =     2.0 # theta pert magnitude
prob.x_c    = 10000.0
prob.z_c    =  2000.0
prob.x_r    =  2000.0
prob.z_r    =  2000.0
prob.T_0    =   300.0
prob.do_moist_bubble = false
prob.T_pert_is_airtemp = false # Perturb theta
