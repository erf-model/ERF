# ------------------  INPUTS TO MAIN PROGRAM  -------------------
stop_time = 10.

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# Smaller for debugging
geometry.prob_extent    =  500.    500.   1000.
amr.n_cell              =   48      48      96

geometry.is_periodic = 1 1 0

zhi.type = "SlipWall"
zlo.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 1.0
erf.fixed_mri_dt_ratio = 4

# DIAGNOSTICS & VERBOSITY
erf.sum_interval    = 1     # timesteps between computing mass
erf.v               = 1     # verbosity in ERF.cpp
amr.v               = 1     # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0     # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk   # root name of checkpoint file
erf.check_int       = 100   # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt   # prefix of plotfile name
erf.plot_int_1      = 1     # number of timesteps between plotfiles
erf.plot_vars_1     = density x_velocity y_velocity z_velocity pressure theta rhoKE pres_hse dens_hse pert_pres pert_dens qv qc

# SOLVER CHOICES
erf.use_gravity = true
erf.molec_diff_type = "None"
erf.les_type = "Smagorinsky"
erf.Cs       = 0.25

erf.moisture_model  = "Kessler"
erf.buoyancy_type   = 1

erf.rayleigh_damp_U = true
erf.rayleigh_damp_V = true
erf.rayleigh_damp_W = true
erf.rayleigh_damp_T = true
    
erf.init_type = "input_sounding"
erf.init_sounding_ideal = true

# PROBLEM PARAMETERS
# these are zeroed because we are using an input_sounding
prob.rho_0 = 0.0
prob.T_0 = 0.0
prob.A_0 = 0.0
prob.U_0 = 0.0
prob.V_0 = 0.0
prob.W_0 = 0.0

prob.pert_ref_height = -1.0
prob.U_0_Pert_Mag = 0.0 # this causes rho < 0 in the first slow step
prob.V_0_Pert_Mag = 0.0
prob.W_0_Pert_Mag = 0.0

prob.pert_deltaU = 0.0
prob.pert_deltaV = 0.0

prob.pert_periods_U = 1
prob.pert_periods_V = 1
