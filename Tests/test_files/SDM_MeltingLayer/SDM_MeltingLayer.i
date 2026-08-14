# ------------------  1D MELTING LAYER (ctest)  --------------------------------
# Monodisperse low-density ice flakes (D = 10 mm, rho_app = 20 kg/m^3) are injected at
# the top of a prescribed, frozen melting layer: a linear T(z) with the 0 C level at the
# top (z = 800 m) warming downward to +4.8 C (lapse 0.6 C/100m, RH 100%). The flakes fall
# under gravity and melt (Seifert-Beheng) as they descend into warmer air; once a flake
# is a mixed ice-water particle its terminal velocity blends from ice-like to water-like
# (Frick et al. 2013). This test therefore exercises melting and the mixed-phase blended
# terminal velocity together. The final plotfile is compared to the gold.
# Deterministic: monodisperse, deterministic placement, fixed RNG seed.

erf.prob_name      = "SDM_Congestus3D_cold"

max_step  = 300
stop_time = 1.0e9
amrex.fpe_trap_invalid = 1
erf.fix_random_seed = 1
fabarray.mfiter_tile_size = 1024 1024 1024

# 8x8 horizontal (domain doubled vs the thin column) so the grid decomposes over 4 ranks
geometry.prob_lo     =  0.0    0.0    0.0
geometry.prob_hi     =  160.0  160.0  800.0
amr.n_cell           =  8      8      80      # dx = dy = 20 m, dz = 10 m
geometry.is_periodic =  1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.fixed_dt       = 0.5
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
erf.plot_int_1   = 300
erf.plot_vars_1  = density temp theta pressure qv qi qrain rel_humidity \
                   super_droplets_moisture_number_density \
                   super_droplets_moisture_species_mass_density_H2O \
                   super_droplets_moisture_species_mass_density_ice
particles.disable_plt = true

# Super-droplet method: prescribed frozen environment, melting only (no coalescence)
super_droplets_moisture.stable_redistribute    = true
super_droplets_moisture.kinematic_mode         = true
super_droplets_moisture.initial_distribution_type = "uniform"
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.diagnostics_interval   = -1
super_droplets_moisture.include_phase_change    = true
super_droplets_moisture.include_coalescence     = false
super_droplets_moisture.advect_with_flow       = true
super_droplets_moisture.advect_with_gravity    = true
super_droplets_moisture.prescribed_advection   = false
super_droplets_moisture.density_scaling        = false
super_droplets_moisture.aerosols               = NaCl

# No initial particles -- the column is populated by injection at the top (the 0 C level)
super_droplets_moisture.num_initializations    = 0

super_droplets_moisture.num_injections = 1
super_droplets_moisture.injection.distribution_type = "uniform"
super_droplets_moisture.injection.multiplicity_type = "constant"
super_droplets_moisture.injection.domain_velocity   = 0.0 0.0 0.0
super_droplets_moisture.injection.particle_box_lo   = 0.0    0.0    790.0
super_droplets_moisture.injection.particle_box_hi   = 160.0  160.0  800.0
super_droplets_moisture.injection.ice_apparent_density = 20.0
super_droplets_moisture.injection.species_distribution_type_ice = "mass_constant"
super_droplets_moisture.injection.species_mean_mass_ice         = 1.047e-05   # D = 10 mm, rho = 20
super_droplets_moisture.injection.species_distribution_type_H2O = "mass_constant"
super_droplets_moisture.injection.species_mean_mass_H2O         = 0.0
super_droplets_moisture.injection.aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.injection.aerosol_mean_mass_NaCl         = 1.0e-19
super_droplets_moisture.injection.rate              = 5.0e1
super_droplets_moisture.injection.particles_per_cell = 1
