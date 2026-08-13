# =============================================================================
#  SDM 2D bubble + fractional injection  (multi-rank: 1 / 2 / 4)
#
#  Based on the ERF regression test SDM_Bubble2D_Adv_wInjection, with a single
#  injection whose rate is lowered so the requested super-droplets per cell per
#  step are sub-unity, exercising fractional injection:
#
#     rate * dt * cell_volume / particles_per_cell
#       = 2e-7 * 0.5 * (100*100*100) / 1 = 0.1   (< 1)
#
#  -> one super-droplet (multiplicity 1) per box cell every 10 steps. The
#  per-injection accumulator is a scalar advanced identically on every rank, so
#  the result is identical for 1, 2, and 4 ranks.
# =============================================================================

erf.prob_name = "Bubble"

max_step  = 50
stop_time = 3600.0

erf.init_type = MoistBaseState
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_extent = 20000.0 400.0 10000.0
amr.n_cell           = 200     4      100        # dx=dy=dz=100 m -> cell volume 1e6 m^3
geometry.is_periodic = 0 1 0
xlo.type = "SlipWall"
xhi.type = "SlipWall"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.fixed_dt = 0.5
erf.fixed_mri_dt_ratio = 4

erf.sum_interval = 1
erf.v            = 1
amr.v            = 1
amr.max_level    = 0
erf.check_int    = -1

erf.plot_file_1 = plt
erf.plot_int_1  = 50
erf.plot_vars_1 = density rhotheta temp pressure \
                  super_droplets_moisture_number_density \
                  super_droplets_moisture_sd_number_density \
                  super_droplets_moisture_aerosol_mass_density_NaCl
particles.disable_plt = true

erf.use_gravity = true
erf.dycore_horiz_adv_type    = "Upwind_3rd"
erf.dycore_vert_adv_type     = "Upwind_3rd"
erf.dryscal_horiz_adv_type   = "Upwind_3rd"
erf.dryscal_vert_adv_type    = "Upwind_3rd"
erf.moistscal_horiz_adv_type = "Upwind_3rd"
erf.moistscal_vert_adv_type  = "Upwind_3rd"

erf.les_type       = "None"
erf.pbl_type       = "None"
erf.moisture_model = "SuperDroplets"
erf.buoyancy_type  = 1
erf.molec_diff_type = "ConstantAlpha"
erf.rho0_trans      = 1.0
erf.alpha_T         = 0.0
erf.alpha_C         = 0.0

# ---- SDM: injection only (no processes), aerosol super-droplets -------------
super_droplets_moisture.stable_redistribute     = true
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.include_phase_change     = false
super_droplets_moisture.include_coalescence      = false
super_droplets_moisture.advect_with_gravity      = false
super_droplets_moisture.distribution_type        = "uniform"
super_droplets_moisture.diagnostics_interval     = -1
super_droplets_moisture.multiplicity_type        = "constant"
super_droplets_moisture.recycle_particles        = false
super_droplets_moisture.aerosols                 = NaCl

super_droplets_moisture.num_initializations        = 0
super_droplets_moisture.initial_number_density     = 0.0
super_droplets_moisture.initial_particles_per_cell = 0

# ---- Single fractional injection (sub-unity per-step multiplicity) ----------
super_droplets_moisture.num_injections = 1
super_droplets_moisture.injection.distribution_type = "uniform"
super_droplets_moisture.injection.domain_velocity   = 18.0 0.0 0.0
super_droplets_moisture.injection.particle_box_lo   =  1000.0   0.0 5000.0
super_droplets_moisture.injection.particle_box_hi   =  1400.0 400.0 5400.0
super_droplets_moisture.injection.aerosol_distribution_type_NaCl = "mass_constant"
super_droplets_moisture.injection.aerosol_mean_mass_NaCl = 1.0e-19
super_droplets_moisture.injection.rate               = 2.0e-7   # sub-unity: 0.1 SD/cell/step
super_droplets_moisture.injection.particles_per_cell = 1

# warm moist bubble (provides the flow field)
prob.x_c = 10000.0
prob.z_c =  2000.0
prob.x_r =  2000.0
prob.z_r =  2000.0
prob.T_0 =   300.0
prob.do_moist_bubble = true
prob.theta_pert  = 2.0
prob.qt_init     = 0.02
prob.eq_pot_temp = 320.0
