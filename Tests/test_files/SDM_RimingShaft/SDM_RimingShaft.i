# ------------------  1D RIMING SHAFT (ctest)  ---------------------------------
# A cold (-10 C), water-saturated, supercooled column is seeded with a uniform cloud
# of monodisperse 10 um droplets (the LWC). Solid ice spheres (D = 2 mm) are injected
# at the top, fall under gravity, and grow by collecting droplets (riming = ice-water
# collection); phase_change is OFF so the only process is collection. A short run
# exercises the riming routine -- ice mass grows as cloud water is swept up.
#
# Collision (and injection) sampling is stochastic; reproducibility relies on
# erf.fix_random_seed + stable_redistribute, so the gold is platform-specific. Gated
# behind ERF_TEST_ENABLE_EXTRA_SDM_TESTS and skipped on GPU.

erf.prob_name                = "SDM_Congestus3D_cold"

max_step  = 400
stop_time = 1.0e9
amrex.fpe_trap_invalid = 1
erf.fix_random_seed = 1
fabarray.mfiter_tile_size = 1024 1024 1024

# 8x8 horizontal (domain doubled vs the thin column) so the grid decomposes over 4 ranks
geometry.prob_lo     =  0.0    0.0    0.0
geometry.prob_hi     =  160.0  160.0  1000.0
amr.n_cell           =  8      8      100     # dx = dy = 20 m, dz = 10 m
geometry.is_periodic =  1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.fixed_dt       = 0.25
erf.init_type      = input_sounding
erf.use_gravity    = true
erf.buoyancy_type  = 2
erf.moisture_model = "SuperDroplets"
erf.les_type       = "None"
erf.pbl_type       = "None"

amr.max_level    = 0
erf.sum_interval = -1
erf.v            = 1
amr.v            = 1
erf.check_int    = -1
erf.plot_file_1  = plt
erf.plot_int_1   = 400
erf.plot_vars_1  = density temp theta pressure qv qc qi rel_humidity \
                   super_droplets_moisture_number_density \
                   super_droplets_moisture_species_mass_density_ice \
                   super_droplets_moisture_species_mass_density_H2O
particles.disable_plt = true

# Super-droplet method: prescribed environment, collection (riming) only
super_droplets_moisture.stable_redistribute    = true
super_droplets_moisture.kinematic_mode         = true
super_droplets_moisture.initial_distribution_type = "uniform"
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.diagnostics_interval   = -1
super_droplets_moisture.include_phase_change    = false
super_droplets_moisture.include_coalescence     = true
super_droplets_moisture.kernel_relative_velocity = "absolute_velocity"
super_droplets_moisture.coalescence_bin_size    = 4 4 1
super_droplets_moisture.advect_with_flow        = true
super_droplets_moisture.advect_with_gravity     = true
super_droplets_moisture.prescribed_advection    = false
super_droplets_moisture.density_scaling         = false
super_droplets_moisture.aerosols                = NaCl

# Seeded supercooled cloud: monodisperse 10 um droplets (LWC 0.7 g/m^3)
super_droplets_moisture.num_initializations    = 1
super_droplets_moisture.multiplicity_type      = "constant"
super_droplets_moisture.initial_number_density = 1.671e8
super_droplets_moisture.initial_particles_per_cell = 10
super_droplets_moisture.initial_species_distribution_type_H2O = "mass_constant"
super_droplets_moisture.initial_species_mean_mass_H2O = 4.189e-12   # 10 um radius droplet
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.initial_aerosol_mean_mass_NaCl = 1.0e-19

# Steady ice flux from the top
super_droplets_moisture.num_injections = 1
super_droplets_moisture.injection.distribution_type = "uniform"
super_droplets_moisture.injection.multiplicity_type = "constant"
super_droplets_moisture.injection.domain_velocity   = 0.0 0.0 0.0
super_droplets_moisture.injection.particle_box_lo   = 0.0  0.0  990.0
super_droplets_moisture.injection.particle_box_hi   = 20.0 20.0 1000.0
super_droplets_moisture.injection.ice_apparent_density = 916.8     # solid ice seed
super_droplets_moisture.injection.species_distribution_type_ice = "mass_constant"
super_droplets_moisture.injection.species_mean_mass_ice         = 3.8418e-06   # D = 2 mm solid
super_droplets_moisture.injection.species_distribution_type_H2O = "mass_constant"
super_droplets_moisture.injection.species_mean_mass_H2O         = 0.0
super_droplets_moisture.injection.aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.injection.aerosol_mean_mass_NaCl         = 1.0e-19
super_droplets_moisture.injection.rate              = 5.0e1
super_droplets_moisture.injection.particles_per_cell = 1
