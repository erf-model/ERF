# ------------------  INPUTS TO MAIN PROGRAM  -------------------
erf.prob_name = "ABL"

max_step = 10

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  40   40  1000
amr.n_cell           =   4    4    51

geometry.is_periodic = 1 1 0

# MOST BOUNDARY (DEFAULT IS ADIABATIC FOR THETA)
zlo.type        = "surface_layer"
erf.surface_layer.flux_type = "custom"
erf.most.ustar  = 0.1
erf.most.tstar  = 0.1  # actually ustar*tstar
erf.most.qstar  = 0.001   # actually ustar*qstar
erf.most.zref   = 10.0

zhi.type        = "SlipWall"

# INITIALIZATION
erf.init_type          = "input_sounding"
erf.input_sounding_file = "input_sounding"

# TIME STEP CONTROL
erf.fixed_mri_dt_ratio = 6
erf.fixed_dt = 10.0
   
# DIAGNOSTICS & VERBOSITY
erf.sum_interval    = 1       # timesteps between computing mass
erf.v               = 1       # verbosity in ERF.cpp
amr.v               = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = -1         # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt       # prefix of plotfile name
erf.plot_int_1      = 10        # number of timesteps between plotfiles
erf.plot_vars_1     = density rhoadv_0 x_velocity y_velocity z_velocity pressure temp theta qv
    
# SOLVER CHOICE
erf.vert_implicit_fac  = 1.0 1.0 0.0  # implicit in first two stages, explicit in last RK stage

erf.moisture_model = SatAdj
  
erf.molec_diff_type = "Constant"
erf.alpha_T         = 10.0 # [m^2/s]
erf.alpha_C         = 10.0
erf.dynamic_viscosity = 10.0

erf.use_gravity  = true
erf.use_coriolis = false

# Turbulence closure
erf.les_type    = "None"
erf.pbl_type    = "None"
