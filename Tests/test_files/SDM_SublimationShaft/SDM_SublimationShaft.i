# ------------------  1D SUBLIMATION SHAFT (ctest)  ----------------------------
# Monodisperse solid-ice spheres (D = 0.35 mm) are seeded through a prescribed,
# frozen, isothermal (T = -15 C), subsaturated-over-ice (RH_i = 82%) column. They fall
# under gravity and lose mass by vapour diffusion (MassChange_SV). A short run exercises
# the sublimation routine; the final plotfile is compared to the gold.
# Deterministic: monodisperse, deterministic placement, fixed RNG seed.

erf.prob_name      = "SDM_Congestus3D_cold"
max_step  = 100
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
erf.plot_int_1   = 100
erf.plot_vars_1  = density temp theta pressure qv qi rel_humidity \
                   super_droplets_moisture_number_density \
                   super_droplets_moisture_species_mass_density_ice \
                   super_droplets_moisture_species_mass_density_H2O
particles.disable_plt = true

# Super-droplet method: kinematic (frozen env), sublimation only, no coalescence
super_droplets_moisture.stable_redistribute    = true
super_droplets_moisture.kinematic_mode         = true
super_droplets_moisture.diagnostics_interval   = -1
super_droplets_moisture.include_phase_change   = true
super_droplets_moisture.include_coalescence    = false
super_droplets_moisture.advect_with_flow       = true
super_droplets_moisture.advect_with_gravity    = true
super_droplets_moisture.prescribed_advection   = false
super_droplets_moisture.density_scaling        = false
super_droplets_moisture.aerosols               = NaCl
super_droplets_moisture.dimensionality         = one_d_z

# Monodisperse initial particles, deterministic placement (1 per cell, whole column)
super_droplets_moisture.initial_distribution_type = "uniform"
super_droplets_moisture.place_randomly_in_cells   = false
super_droplets_moisture.multiplicity_type         = "constant"
super_droplets_moisture.ice_apparent_density      = 916.8        # solid ice sphere
super_droplets_moisture.initial_species_distribution_type_ice = "mass_constant"
super_droplets_moisture.initial_species_mean_mass_ice         = 2.058e-08   # D=0.35 mm solid
super_droplets_moisture.initial_species_distribution_type_H2O = "mass_constant"
super_droplets_moisture.initial_species_mean_mass_H2O         = 0.0
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.initial_aerosol_mean_mass_NaCl         = 1.0e-19
super_droplets_moisture.initial_number_density     = 1.0e6
super_droplets_moisture.initial_particles_per_cell = 1
