# Sampled multiplicities of super-droplets with more than one varying component, at t = 0.
# Mode 0: a dust core and its sulfate coating, two log-normal aerosols. Mode 1: an exponential
# water mass on a constant aerosol. The multiplicity-weighted moments of each component must
# match its distribution, and each mode its number density (check_sampled_multiplicity.py).
max_step  = 0
stop_time = 1.0

amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo     = 0. 0. 0.
geometry.prob_hi     = 400. 400. 400.
amr.n_cell           = 8 8 8
geometry.is_periodic = 1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.prob_name = "Bubble"
prob.T_pert   = 0.0
erf.init_type = "uniform"
erf.use_gravity = false

erf.fixed_dt         = 1.0
erf.substepping_type = "None"

erf.sum_interval = -1
erf.v            = 1
amr.max_level    = 0
erf.check_int    = -1
erf.plot_file_1  = plt
erf.plot_int_1   = 1
erf.plot_vars_1  = density super_droplets_moisture_number_density
particles.disable_plt = false

erf.moisture_model = "SuperDroplets"
erf.random_seed    = 7
super_droplets_moisture.initial_distribution_type    = "uniform"
super_droplets_moisture.diagnostics_interval         = -1
super_droplets_moisture.include_phase_change         = false
super_droplets_moisture.include_advection            = false
super_droplets_moisture.include_coalescence          = false
super_droplets_moisture.multiplicity_type            = "sampled"
super_droplets_moisture.aerosols                     = soil NH4HSO4
super_droplets_moisture.num_initializations          = 2

super_droplets_moisture.0.initial_aerosol_distribution_type_soil    = "radius_log_normal"
super_droplets_moisture.0.initial_aerosol_mean_radius_soil          = 2.5e-6
super_droplets_moisture.0.initial_aerosol_geomstd_radius_soil       = 1.70
super_droplets_moisture.0.initial_aerosol_min_radius_soil           = 5.0e-7
super_droplets_moisture.0.initial_aerosol_max_radius_soil           = 1.5e-5
super_droplets_moisture.0.initial_aerosol_distribution_type_NH4HSO4 = "radius_log_normal"
super_droplets_moisture.0.initial_aerosol_mean_radius_NH4HSO4       = 80.0e-9
super_droplets_moisture.0.initial_aerosol_geomstd_radius_NH4HSO4    = 1.50
super_droplets_moisture.0.initial_aerosol_min_radius_NH4HSO4        = 1.0e-9
super_droplets_moisture.0.initial_aerosol_max_radius_NH4HSO4        = 1.0e-6
super_droplets_moisture.0.initial_number_density                    = 1.0e6
super_droplets_moisture.0.initial_particles_per_cell                = 128

super_droplets_moisture.1.initial_aerosol_distribution_type_NH4HSO4 = "mass_constant"
super_droplets_moisture.1.initial_aerosol_mean_mass_NH4HSO4         = 1.0e-19
super_droplets_moisture.1.initial_species_distribution_type_H2O     = "mass_exponential"
super_droplets_moisture.1.initial_species_mean_mass_H2O             = 1.0e-10
super_droplets_moisture.1.initial_species_min_mass_H2O              = 1.0e-14
super_droplets_moisture.1.initial_species_max_mass_H2O              = 5.0e-10
super_droplets_moisture.1.initial_number_density                    = 1.0e8
super_droplets_moisture.1.initial_particles_per_cell                = 128
