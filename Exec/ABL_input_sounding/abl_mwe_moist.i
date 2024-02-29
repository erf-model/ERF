# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#stop_time = 86400.0
stop_time = 115200.0

stop_time = 600.

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# Full size
geometry.prob_extent    = 1000.   1000.   1000.
amr.n_cell              =   96      96      96

# Smaller for debugging
geometry.prob_extent    =  500.    500.   1000.
amr.n_cell              =   48      48      96

geometry.is_periodic = 1 1 0

zhi.type = "SlipWall"
zhi.theta_grad   = 0.003    # for case N02

# MOST BOUNDARY (DEFAULT IS ADIABATIC FOR THETA)
#zlo.type                = "Most"
#erf.most.z0             = 0.01
#erf.most.zref           = 5.21 # >=dz/2
#erf.most.surf_temp_flux = 0.0 # [K-m/s] for case N02

zlo.type = "SlipWall"

erf.use_rayleigh_damping = true

# TIME STEP CONTROL
erf.cfl             = 0.5
erf.no_substepping  = 0

# DIAGNOSTICS & VERBOSITY
erf.sum_interval    = 1     # timesteps between computing mass
erf.v               = 1     # verbosity in ERF.cpp
amr.v               = 1     # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       =  100     # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt       # prefix of plotfile name
erf.plot_int_1      = 1    # number of timesteps between plotfiles
erf.plot_vars_1     = density x_velocity y_velocity z_velocity pressure theta rhoKE pres_hse dens_hse pert_pres pert_dens qv

# SOLVER CHOICES
erf.dycore_horiz_adv_type    = "Centered_2nd"
erf.dycore_vert_adv_type     = "Centered_2nd"
erf.dryscal_horiz_adv_type   = "Centered_2nd"
erf.dryscal_vert_adv_type    = "Centered_2nd"

erf.use_gravity = true
erf.use_gravity = false

#erf.use_coriolis = true
#erf.latitude = 90.0
#erf.rotational_time_period = 125663.706143592  # to get f=1e-4 1/s

#erf.abl_driver_type = "GeostrophicWind"
#erf.abl_geo_wind = 10. 0. 0.

erf.molec_diff_type = "None"
erf.les_type = "Deardorff"
erf.Ck       = 0.1
erf.sigma_k  = 1.0
erf.Ce       = 0.1
erf.RhoKE_0  = 0.1

erf.moisture_model  = "Kessler"
erf.init_type = "input_sounding"
erf.input_sounding_file = "input_sounding_moist"
erf.init_sounding_ideal = true

# PROBLEM PARAMETERS
# these are zeroed because we are using an input_sounding
prob.rho_0 = 0.0
prob.T_0 = 0.0
prob.A_0 = 0.0
prob.U_0 = 0.0
prob.V_0 = 0.0
prob.W_0 = 0.0

prob.pert_ref_height = 25.0
prob.U_0_Pert_Mag = 0.0 # this causes rho < 0 in the first slow step
prob.V_0_Pert_Mag = 0.0
prob.W_0_Pert_Mag = 0.0

prob.pert_deltaU = 1.0
prob.pert_deltaV = 1.0

prob.pert_periods_U = 1
prob.pert_periods_V = 1
