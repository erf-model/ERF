# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression motivation:
# SuperDroplets stores water qv, qc, and qr in RhoQ1-RhoQ3 after readInputs().
# This test prevents the constructor's temporary nonmoist sentinel from being
# misread as a production layout that omits RhoQ3, qrain, qt, or qp.
erf.prob_name = "Bubble"
max_step = 0
stop_time = 5.0
erf.init_type = Uniform

amrex.fpe_trap_invalid = 0
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo = 0. 0. 0.
geometry.prob_hi = 0.8 0.8 0.8
amr.n_cell = 8 8 8
geometry.is_periodic = 1 1 1

erf.fixed_dt = 0.0002
erf.substepping_type = "None"
erf.sum_interval = 1
erf.v = 1
amr.v = 1
amr.max_level = 0

erf.check_file = chk
erf.check_int = -1

erf.plot_file_1 = plt
erf.plot_int_1 = 1
erf.plot_vars_1 = density theta qv qc qrain qt qn qp

erf.use_gravity = false
erf.les_type = "None"
erf.pbl_type = "None"
erf.moisture_model = "SuperDroplets"
erf.buoyancy_type = 1

super_droplets_moisture.stable_redistribute = true
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.distribution_type = "uniform"
super_droplets_moisture.diagnostics_interval = 1
super_droplets_moisture.advect_with_flow = false
super_droplets_moisture.advect_with_gravity = false
super_droplets_moisture.include_coalescence = false
super_droplets_moisture.aerosols = NaCl
super_droplets_moisture.multiplicity_type = "constant"
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.initial_aerosol_mean_mass_NaCl = 1.0e-19
super_droplets_moisture.initial_species_distribution_type_H2O = "mass_constant"
super_droplets_moisture.initial_species_mean_mass_H2O = 4.1887902e-6
super_droplets_moisture.initial_number_density = 1.0e7
super_droplets_moisture.initial_particles_per_cell = 1

prob.U_0 = 0.0
prob.x_c = 2.0
prob.z_c = 2.0
prob.x_r = 1e-99
prob.z_r = 1e-99
prob.do_moist_bubble = false
prob.theta_pert = 0.0
prob.qt_init = 0.02
prob.eq_pot_temp = 300.0
