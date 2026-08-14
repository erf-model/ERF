# ------------------  0D ICE-AGGREGATION BOX (ctest)  --------------------------
# A well-mixed box of polydisperse ice super-droplets (Brdar-Seifert 2018, Sec 2.4.1)
# self-aggregate by ice-ice collection (one coalescence bin = the whole box).
# Aggregation reduces the number density while conserving mass. A short run exercises
# the aggregation routine; the final plotfile is compared to the gold.
#
# Collision (and mass_exponential) sampling is stochastic; reproducibility relies on
# erf.fix_random_seed + stable_redistribute, so the gold is platform-specific. Gated
# behind ERF_TEST_ENABLE_EXTRA_SDM_TESTS and skipped on GPU.

max_step  = 4500
stop_time = 1.0e9
amrex.fpe_trap_invalid = 1
erf.fix_random_seed = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo     =  0.   0.   0.
geometry.prob_hi     =  200. 200. 100.   # 8x8x4 -> 4 MPI ranks; dx=dy=dz=25 m
amr.n_cell           =  8    8    4
geometry.is_periodic =  1 1 1

erf.prob_name        = "Bubble"
erf.fixed_dt         = 1.0
erf.substepping_type = "None"
erf.init_type        = "uniform"
erf.use_gravity      = false
erf.moisture_model   = "SuperDroplets"

amr.max_level    = 0
erf.sum_interval = -1
erf.v            = 1
amr.v            = 1
erf.check_int    = -1
erf.plot_file_1  = plt
erf.plot_int_1   = 4500
erf.plot_vars_1  = density \
                   super_droplets_moisture_number_density \
                   super_droplets_moisture_mass_density \
                   super_droplets_moisture_species_mass_density_ice
particles.disable_plt = true

# Super-droplet method: ice self-aggregation only
super_droplets_moisture.stable_redistribute          = true
super_droplets_moisture.initial_distribution_type    = "uniform"
super_droplets_moisture.place_randomly_in_cells       = false
super_droplets_moisture.diagnostics_interval          = -1
super_droplets_moisture.include_phase_change          = false
super_droplets_moisture.include_coalescence           = true
super_droplets_moisture.kernel_relative_velocity = "absolute_velocity"
super_droplets_moisture.include_brownian_coalescence  = false
super_droplets_moisture.coalescence_bin_size          = 4 4 4
super_droplets_moisture.advect_with_flow              = false
super_droplets_moisture.advect_with_gravity           = false
super_droplets_moisture.multiplicity_type             = "constant"
super_droplets_moisture.aerosols                      = NaCl
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.initial_aerosol_mean_mass_NaCl = 1.0e-22
super_droplets_moisture.initial_species_distribution_type_ice = "mass_exponential"
super_droplets_moisture.initial_species_mean_mass_ice = 1.19047619e-9
super_droplets_moisture.ice_apparent_density          = 100.0
super_droplets_moisture.initial_number_density        = 8.4e5
super_droplets_moisture.initial_particles_per_cell    = 64
