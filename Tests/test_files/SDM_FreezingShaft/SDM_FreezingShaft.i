# ------------------  1D FREEZING SHAFT (ctest)  -------------------------------
# Monodisperse supercooled droplets (D = 20 um), each carrying a monodisperse
# insoluble INP ("soil"), are seeded through a prescribed, frozen column that cools
# with height (-8 C at the base to -40 C aloft, water-saturated). Each droplet is
# assigned a freezing temperature by inverse-CDF sampling of the Niemand (2012) INAS
# curve and freezes (MassChange_SL) where the local temperature is below it, so the
# frozen fraction varies with height along the INAS curve. The singular INAS freezing is
# instantaneous, so a single step captures it; stopping at step 1 also isolates the
# freezing from the deposition growth that follows on the frozen particles in the
# water-saturated (ice-supersaturated) column. The final plotfile is compared to the gold.
#
# Stochastic Tfz sampling is made reproducible by erf.fix_random_seed; the gold is
# platform-specific, so this test is gated behind ERF_TEST_ENABLE_EXTRA_SDM_TESTS.

erf.prob_name      = "SDM_Congestus3D_cold"
max_step  = 1
stop_time = 1.0e9
amrex.fpe_trap_invalid = 1
erf.fix_random_seed = 1
fabarray.mfiter_tile_size = 1024 1024 1024

# 8x8 horizontal (domain doubled vs the thin column) so the grid decomposes over 4 ranks
geometry.prob_lo     =  0.0    0.0    0.0
geometry.prob_hi     =  160.0  160.0  4000.0
amr.n_cell           =  8      8      100     # dx = dy = 20 m, dz = 40 m
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
erf.plot_int_1   = 1
erf.plot_vars_1  = density temp theta pressure qv qc rel_humidity \
                   super_droplets_moisture_number_density \
                   super_droplets_moisture_species_mass_density_ice \
                   super_droplets_moisture_species_mass_density_H2O
particles.disable_plt = true

# Super-droplet method: kinematic (frozen env), immersion freezing only
super_droplets_moisture.stable_redistribute    = true
super_droplets_moisture.kinematic_mode         = true
super_droplets_moisture.diagnostics_interval   = -1
super_droplets_moisture.include_phase_change   = true
super_droplets_moisture.include_coalescence    = false
super_droplets_moisture.advect_with_flow       = true
super_droplets_moisture.advect_with_gravity    = false
super_droplets_moisture.prescribed_advection   = false
super_droplets_moisture.density_scaling        = false
super_droplets_moisture.aerosols               = soil

# Monodisperse initial droplets + INP, deterministic placement (1 per cell)
super_droplets_moisture.initial_distribution_type = "uniform"
super_droplets_moisture.place_randomly_in_cells   = false
super_droplets_moisture.multiplicity_type         = "constant"
super_droplets_moisture.initial_species_distribution_type_H2O = "mass_constant"
super_droplets_moisture.initial_species_mean_mass_H2O         = 4.189e-12   # D=20 um droplet
super_droplets_moisture.initial_species_distribution_type_ice = "mass_constant"
super_droplets_moisture.initial_species_mean_mass_ice         = 0.0
super_droplets_moisture.initial_aerosol_distribution_type_soil = "mass_constant"
super_droplets_moisture.initial_aerosol_mean_mass_soil         = 4.0e-13   # INP -> A ~ 2.3e-10 m^2
super_droplets_moisture.initial_number_density     = 1.0e6
super_droplets_moisture.initial_particles_per_cell = 1
