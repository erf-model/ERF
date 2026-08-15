#------------------  INPUTS TO MAIN PROGRAM  -------------------
erf.prob_name = "ABL"

max_step = 20

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  2    2   1
amr.n_cell           =  64  64  32

geometry.is_periodic = 0 0 0

# MOST on all sides
zlo.type = "surface_layer"
zhi.type = "surface_layer"
xlo.type = "surface_layer"
xhi.type = "surface_layer"
ylo.type = "surface_layer"
yhi.type = "surface_layer"

erf.xlo.surface_layer.flux_type = MOENG 
erf.ylo.surface_layer.flux_type = MOENG
erf.zlo.surface_layer.flux_type = MOENG
erf.xhi.surface_layer.flux_type = MOENG
erf.yhi.surface_layer.flux_type = MOENG
erf.zhi.surface_layer.flux_type = MOENG

erf.zlo.most.surf_temp = 299.0
erf.zhi.most.surf_temp = 280.0
erf.xlo.most.surf_temp = 289.5
erf.xhi.most.surf_temp = 289.5
erf.ylo.most.surf_temp = 289.5
erf.yhi.most.surf_temp = 289.5

erf.zlo.most.surf_moist = 0.02143 
erf.zhi.most.surf_moist = 0.00623 
erf.xlo.most.surf_moist = 0.00765 # 65% saturation 
erf.xhi.most.surf_moist = 0.00765 
erf.ylo.most.surf_moist = 0.00765 
erf.yhi.most.surf_moist = 0.00765

erf.zlo.most.z0 = 0.000035
erf.zhi.most.z0 = 0.000035
erf.xlo.most.z0 = 0.000035
erf.xhi.most.z0 = 0.000035
erf.ylo.most.z0 = 0.000035
erf.yhi.most.z0 = 0.000035

erf.xlo.most.zref = 0.03125
erf.ylo.most.zref = 0.03125
erf.zlo.most.zref = 0.03125
erf.xhi.most.zref = 0.03125
erf.yhi.most.zref = 0.03125
erf.zhi.most.zref = 0.03125

# Regional MOST at point
erf.zlo.most.average_policy = 1
erf.zlo.most.radius = 0
erf.zhi.most.average_policy = 1 
erf.zhi.most.radius = 0
erf.xlo.most.average_policy = 1
erf.xlo.most.radius = 0
erf.xhi.most.average_policy = 1
erf.xhi.most.radius = 0
erf.ylo.most.average_policy = 1
erf.ylo.most.radius = 0
erf.yhi.most.average_policy = 1
erf.yhi.most.radius = 0

# Planar MOST
#erf.zlo.most.average_policy = 0
#erf.zhi.most.average_policy = 0
#erf.xlo.most.average_policy = 0
#erf.xhi.most.average_policy = 0
#erf.ylo.most.average_policy = 0
#erf.yhi.most.average_policy = 0

erf.vert_implicit = false
erf.terrain_type = StaticFittedMesh

# TIME STEP CONTROL
erf.anelastic = 1
#erf.use_fft = 1 
erf.fixed_dt           = 0.025  # fixed time step depending on grid resolution
erf.fixed_mri_dt_ratio = 4

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 0       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk
erf.check_int       = 20 number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt
erf.plot_int_1      = 10

erf.plot_vars_1     = x_velocity y_velocity z_velocity temp theta qv qsat qc qt pressure density 
#erf.plot_vars_1     = x_velocity y_velocity z_velocity temp theta qv qsat qc qt pressure density Tau11 Tau12 Tau13 Tau21 Tau22 Tau23 Tau31 Tau32 Tau33 hfx1 hfx2 hfx3 q1fx1 q1fx2 q1fx3 

# SOLVER CHOICE
erf.alpha_T = 0
erf.alpha_C = 0
erf.use_gravity = true 

erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.17
erf.buoyancy_type   = 1

# Initialization
erf.init_type           = "input_sounding"
erf.init_sounding_ideal = true
erf.input_sounding_file = "input_sounding"

# moisture model
erf.moisture_model = "MoistNoCondensation"

#PROBLEM PARAMETERS
prob.rho_0 = 1.2

prob.U_0 = 0.0
prob.V_0 = 0.0
prob.W_0 = 0.0
prob.T_0 = 289.5

# Higher values of perturbations lead to instability
# Instability seems to be coming from BC
prob.T_0_Pert_Mag = 0.
prob.pert_rhotheta = false
prob.U_0_Pert_Mag = 0.05
