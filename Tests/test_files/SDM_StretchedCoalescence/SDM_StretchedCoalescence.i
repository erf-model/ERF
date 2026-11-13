# Golovin coalescence on a stretched mesh: four 50 m layers under four 150 m layers
# (computational spacing 100 m), no advection or phase change. Every layer must start
# at the specified number density and lose droplets at the Golovin rate b*phi of its
# own liquid volume fraction phi (check_stretched_coalescence.py).
max_step  = 300
stop_time = 300.0

amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo     = 0. 0. 0.
geometry.prob_hi     = 200. 200. 800.
amr.n_cell           = 4 4 8
geometry.is_periodic = 1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"
erf.terrain_z_levels = 0.0 50.0 100.0 150.0 200.0 350.0 500.0 650.0 800.0

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
erf.plot_int_1   = 300
erf.plot_vars_1  = density super_droplets_moisture_number_density
particles.disable_plt = false

erf.moisture_model = "SuperDroplets"
erf.random_seed    = 7
super_droplets_moisture.initial_distribution_type    = "uniform"
super_droplets_moisture.diagnostics_interval         = -1
super_droplets_moisture.include_phase_change         = false
super_droplets_moisture.include_advection            = false
super_droplets_moisture.include_coalescence          = true
super_droplets_moisture.include_cold_processes       = false
super_droplets_moisture.coalescence_kernel           = "golovin"
super_droplets_moisture.include_brownian_coalescence = false
super_droplets_moisture.coalescence_bin_size         = 1 1 1
super_droplets_moisture.aerosols                     = NaCl
super_droplets_moisture.multiplicity_type            = "constant"
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.initial_aerosol_mean_mass_NaCl         = 1.0e-22
super_droplets_moisture.initial_species_distribution_type_H2O  = "mass_exponential"
super_droplets_moisture.initial_species_mean_mass_H2O          = 1.1920972798965587e-10
super_droplets_moisture.initial_number_density                 = 8388577.147720466
super_droplets_moisture.initial_particles_per_cell             = 128
