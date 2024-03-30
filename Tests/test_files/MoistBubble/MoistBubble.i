# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step  = 10
stop_time = 3600.0

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 20000.0 400.0  10000.0
amr.n_cell           = 200     4      100
geometry.is_periodic = 0 1 0
xlo.type = "SlipWall"
xhi.type = "SlipWall"    
zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 0.5
erf.fixed_mri_dt_ratio = 4
#erf.no_substepping = 1
#erf.fixed_dt = 0.1

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = 100       # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 100        # number of timesteps between plotfiles
erf.plot_vars_1     = density rhotheta rhoQ1 rhoQ2 rhoadv_0 x_velocity y_velocity z_velocity pressure theta scalar temp pres_hse dens_hse pert_pres pert_dens eq_pot_temp qt qv qc 

# SOLVER CHOICES
erf.use_gravity          = true
erf.use_coriolis         = false
erf.use_rayleigh_damping = false
    
erf.dycore_horiz_adv_type    = "Upwind_3rd"
erf.dycore_vert_adv_type     = "Upwind_3rd"
erf.dryscal_horiz_adv_type   = "Upwind_3rd"
erf.dryscal_vert_adv_type    = "Upwind_3rd"
erf.moistscal_horiz_adv_type = "Upwind_3rd"
erf.moistscal_vert_adv_type  = "Upwind_3rd"       

# PHYSICS OPTIONS
erf.les_type        = "None"
erf.pbl_type        = "None"
erf.moisture_model  = "Kessler_NoRain"
erf.buoyancy_type   = 1
erf.use_moist_background = true

erf.molec_diff_type  = "ConstantAlpha"
erf.rho0_trans       = 1.0 # [kg/m^3], used to convert input diffusivities
erf.dynamicViscosity = 0.0 # [kg/(m-s)] ==> nu = 75.0 m^2/s
erf.alpha_T          = 0.0 # [m^2/s]
erf.alpha_C          = 0.0

# INITIAL CONDITIONS
#erf.init_type = "input_sounding"
#erf.input_sounding_file = "BF02_moist_sounding"
#erf.init_sounding_ideal = true

# PROBLEM PARAMETERS (optional)
# warm bubble input
prob.x_c    = 10000.0
prob.z_c    =  2000.0
prob.x_r    =  2000.0
prob.z_r    =  2000.0
prob.T_0    =   300.0

prob.do_moist_bubble = true
prob.theta_pert  = 2.0
prob.qt_init     = 0.02
prob.eq_pot_temp = 320.0
