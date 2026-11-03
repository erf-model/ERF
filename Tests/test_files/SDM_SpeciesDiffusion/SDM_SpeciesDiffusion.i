# A vapor slab in a sheared, constant-density column. Appended condensate
# belongs to particles; only the appended vapors receive host diffusion.
erf.prob_name = "MultiSpecies Bubble"
erf.init_type = input_sounding
erf.input_sounding_file = input_sounding
erf.sounding_type = ConstantDensity
max_step = 4
erf.fixed_dt = 0.02
erf.fixed_mri_dt_ratio = 16
erf.substepping_type = Implicit
erf.vert_implicit_fac = 1
erf.implicit_moisture_diffusion = true
erf.use_gravity = false
erf.transport_scalar = false
amrex.fpe_trap_invalid = 1
erf.fix_random_seed = 1

geometry.prob_extent = 8 4 16
amr.n_cell = 8 4 16
amr.max_level = 0
amr.max_grid_size = 32
amr.max_grid_size_z = 1024
geometry.is_periodic = 1 1 0
zlo.type = SlipWall
zhi.type = SlipWall
erf.coupling_type = TwoWay

erf.v = 1
erf.sum_interval = 1
erf.plot_file_1 = plt
erf.plot_int_1 = 1
erf.plot_vars_1 = density rhotheta rhoQ7 rhoQ8 rhoQ9 rhoQ10 \
                  qv_water qc_water qv_agua qc_agua \
                  super_droplets_moisture_species_mass_density_water \
                  super_droplets_moisture_species_mass_density_agua Kmv
erf.check_int = -1

erf.moisture_model = SuperDroplets
erf.les_type = Smagorinsky
erf.Cs = 12.0
erf.molec_diff_type = None
erf.dycore_horiz_adv_type = Upwind_3rd
erf.dycore_vert_adv_type = Upwind_3rd
erf.moistscal_horiz_adv_type = Upwind_3rd
erf.moistscal_vert_adv_type = Upwind_3rd

prob.do_moist_bubble = false
prob.T_pert = 0.0
prob.x_r = 0
prob.y_r = 0
prob.z_c = 8
prob.z_r = 4
prob.qv_init_water = 0.0025
prob.qv_init_agua = 0.005

super_droplets_moisture.species = water agua
super_droplets_moisture.aerosols = NaCl
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = mass_constant
super_droplets_moisture.initial_aerosol_mean_mass_NaCl = 1.0e-19
super_droplets_moisture.include_cold_processes = false
super_droplets_moisture.include_phase_change = false
super_droplets_moisture.include_advection = false
super_droplets_moisture.include_coalescence = false
super_droplets_moisture.stable_redistribute = true
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.distribution_type = uniform
super_droplets_moisture.initial_number_density = 1.0e7
super_droplets_moisture.initial_particles_per_cell = 1
super_droplets_moisture.multiplicity_type = constant
super_droplets_moisture.initial_species_mean_mass_H2O = 1.0e-12
super_droplets_moisture.initial_species_mean_mass_water = 1.0e-12
super_droplets_moisture.initial_species_mean_mass_agua = 2.0e-12
super_droplets_moisture.diagnostics_interval = 0
