# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# The moist flat.i ABL of ImmersedForm7MoistEquivalence over a heated, evaporating surface: a kinematic
# heat flux of 0.02 K m/s and a moisture flux of 1e-4 kg/kg m/s (about 24 and 300 W/m2), with the moist
# TKE buoyancy of the k-eqn closure (erf.rans_moist_buoyancy) and its length scale from theta_v. The
# reference for plane.i, the same ABL over an immersed plane on a cell face with the fraction-stress
# wall law and the same fluxes (erf.if_surf_temp_flux, erf.if_surf_moist_flux): the wall moisture
# flux, the Obukhov length from theta_v* and the moist buoyancy on the wall face must give the wind,
# theta, TKE and qv of flat ground. Tests/RunImmersedFlatEquivalence.cmake runs both decks and compares
# the horizontally averaged profiles with Tests/ImmersedFlatEquivalenceCheck.cpp.
#
erf.prob_name = "ABL"
stop_time = 1.0e6
max_step  = 100
amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo     = 0.0 0.0 0.0
geometry.prob_hi     = 160.0 160.0 1000.0
amr.n_cell           = 8 8 50
amr.max_grid_size    = 256
amr.blocking_factor  = 2
amr.max_level        = 0
geometry.is_periodic = 1 1 0

zlo.type    = "surface_layer"
erf.most.z0 = 0.1
erf.most.surf_temp_flux  = 0.02      # kinematic heat flux [K m/s]
erf.most.surf_moist_flux = 1.0e-4    # kinematic moisture flux [kg/kg m/s]
erf.dirichlet_k = true    # wall k from u* (AL01 Eq. 16)
zhi.type       = "SlipWall"
zhi.theta_grad = 0.0

# anelastic, implicit vertical diffusion (solved in the first MidPoint stage)
erf.compressibility_type = Anelastic
erf.timestepping_type = MidPoint
erf.vert_implicit  = true
erf.use_fft        = true    # FFT Poisson solve (flat, periodic in x and y)
erf.fixed_dt       = 1.0

erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.sounding_type       = "Ideal"

erf.sum_interval    = 60
erf.v               = 1
amr.v               = 0

erf.check_file  = chk
erf.check_int   = -1
erf.plot_file_1 = plt_flat
erf.plot_int_1  = 100
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta KE Kmv walldist qv

erf.dycore_horiz_adv_type  = "Upwind_5th"
erf.dycore_vert_adv_type   = "Upwind_3rd"
erf.dryscal_horiz_adv_type = "Upwind_5th"
erf.dryscal_vert_adv_type  = "Upwind_3rd"

erf.use_gravity  = true
erf.use_coriolis = true
erf.coriolis_3d  = false
erf.latitude     = 90.0
erf.rotational_time_period = 125663.7061435917   # f = 1e-4 1/s
erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind    = 10.0 0.0 0.0

# top sponge toward the sounding (10 m/s), as the kynema-sgf deck; dampcoef * dt = 0.05
erf.rayleigh_damp_U   = true
erf.rayleigh_damp_V   = true
erf.rayleigh_damp_W   = true
erf.rayleigh_dampcoef = 0.05
erf.rayleigh_zdamp    = 200.0

erf.moisture_model  = "MoistNoCondensation"   # qv and qc transported, no condensation
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.rans_type       = "kEqn"
erf.thermal_stratification = thetav   # k-eqn length scale from the virtual potential temperature
erf.rans_moist_buoyancy    = true     # moist TKE buoyancy, g/theta_v0 [(1 + epsv qv) <w'theta'> + epsv theta <w'qv'>]
erf.theta_ref       = 300.0

# same initial k in every case: 0.84 m2/s2, linear decay to zero at 900 m above z = 0
prob.KE_0            = 0.84
prob.KE_decay_height = 900.0
prob.KE_decay_order  = 1
prob.pert_deltaU  = 0.0
prob.pert_deltaV  = 0.0
prob.T_0_Pert_Mag = 0.0
prob.U_0_Pert_Mag = 0.0
prob.V_0_Pert_Mag = 0.0
prob.W_0_Pert_Mag = 0.0
