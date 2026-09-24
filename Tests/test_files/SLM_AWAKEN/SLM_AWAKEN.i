# ------------------  INPUTS TO MAIN PROGRAM  -------------------
stop_time = 169200
max_step = 25200

amrex.fpe_trap_invalid = 0
fabarray.mfiter_tile_size = 1024 1024 1024

# Set to 1/4 of 40GB (10GB) or 1/3 (13GB) or even smaller
amrex.the_arena_init_size = 5368709120 #10737418240 #can help with memory issues when using radiation

erf.check_for_nans = 2 # for debugging

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 400000.        400000.       17500.
amr.n_cell           = 400            400           47

geometry.is_periodic = 0 0 0

# MOST BOUNDARY (DEFAULT IS ADIABATIC FOR THETA)
xlo.type = "Outflow"
xhi.type = "Outflow"
ylo.type = "Outflow"
yhi.type = "Outflow"

zlo.type = "surface_layer"
erf.most.z0                = 0.001  # for water body features when using SLM
erf.most.surf_temp_flux    = 0.0
erf.most.zref              = 15.0 # half first dz
erf.most.average_policy    = 1
erf.most.radius            = 0

zhi.type        = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 1
erf.substepping_cfl = 0.5
erf.substepping_type = "Implicit"

# DIAGNOSTICS & VERBOSITY
erf.sum_interval    = 1       # timesteps between computing mass
erf.v               = 1       # verbosity in ERF.cpp
amr.v               = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = -1 #3600         # number of timesteps between checkpoints
erf.restart_type    = "native"

# PLOTFILES
erf.plot_file_1     = plt       # prefix of plotfile name
erf.plot_int_1      = 600      # number of timesteps between plotfiles
erf.plot_vars_1     = detJ Kmv x_velocity y_velocity z_velocity temperature theta z_phys qt qn qp qc qi qv

erf.plot2d_file_1  = sfc
erf.plot2d_int_1   = 600
erf.plot2d_vars_1  = t_surf u_star

# SOLVER CHOICE
erf.dycore_horiz_adv_type    = Upwind_3rd
erf.dycore_vert_adv_type     = Upwind_3rd
erf.dryscal_horiz_adv_type   = Upwind_3rd_SL
erf.dryscal_vert_adv_type    = Upwind_3rd_SL
erf.moistscal_horiz_adv_type = Upwind_3rd_SL
erf.moistscal_vert_adv_type  = Upwind_3rd_SL
erf.use_gravity = true

# Turbulence closure
erf.pbl_type         = MRF
prob.KE_0            = 0.4
prob.KE_decay_height = 500.
prob.KE_decay_order  = 3

erf.les_type        = "Smagorinsky2D" #"Smagorinsky"
erf.Cs              = 0.1
#erf.les_type = "Smagorinsky"
#erf.Cs = 0.1
erf.use_smag_stratification = false
erf.mix_isotropic = false

erf.implicit_before_substep     = true
erf.implicit_thermal_diffusion  = true
erf.implicit_moisture_diffusion = true
erf.implicit_momentum_diffusion = true
erf.vert_implicit_fac           = 1. 1. 0.

erf.use_coriolis      = true
erf.variable_coriolis = true
erf.coriolis_3d       = true

# Terrain
erf.terrain_type = "StaticFittedMesh"
erf.terrain_smoothing = 2

# Rayleigh
erf.rayleigh_damping_type = FastImplicit
erf.rayleigh_damp_W   = true
erf.rayleigh_zdamp    = 5000.0
erf.rayleigh_dampcoef = 0.2

# Moisture
erf.moisture_model = "MoistNoCondensation"

# INITIALIZATION WITH ATM DATA
erf.real_width     = 5
erf.real_set_width = 1
erf.init_type      = "WRFInput"
erf.nc_init_file_0 = "./wrfinput_d01"
erf.nc_bdy_file    = "./wrfbdy_d01"

# SLM
erf.land_surface_model = "SLM"
erf.plot_lsm = true           # whether to plot 2D/3D SLM fields
slm.nsoil = 4                 # number of soil layers, must match wrfinput file definition
slm.soil_dz = 0.1 0.3 0.6 1.0 # SLM layer layout; WRFInput DZS supplies per-cell thickness
slm.tabs_s = 300.0            # fallback surface temperature for non-land cells (K)
slm.t00 = 300.0               # constant offset for sstxy

slm.z0_soil = 0.0387  # baresoil roughness length
slm.mws_mx0 = 50.0    # maximum puddle water storage (mm)
slm.Rc_max  = 5000.0  # maximum stomatal resistance
slm.T_opt   = 298.0   # optimum temperature for transpiration
slm.zref    = 15.0    # height of reference level (m)

slm.use_parameter_file = true                # load SLM soil/vegetation parameters from a Noah-MP .TBL file
slm.parameter_file     = "NoahmpTable.TBL"   # TBL file containing parameters for SLM
slm.veg_dataset        = "modis"             # Vegetation parameter dataset to use from the parameter_file
slm.soil_dataset       = "stas"              # Soil parameter dataset to use from the parameter_file
slm.use_param_tbl      = false               # enable WRF/Noah-style LAI and vegetation-fraction handling

# Radiation setting
erf.o3vmr                  = 4.825e-08
erf.radiation_model        = "RRTMGP" # or "None" to disable
erf.nswbands               = 14
erf.nlwbands               = 16
erf.nswgpts                = 112
erf.nlwgpts                = 128
erf.rrtmgp_file_path       = ./
erf.rrtmgp_coeffs_sw       = rrtmgp-gas-sw-g112.nc
erf.rrtmgp_coeffs_lw       = rrtmgp-gas-lw-g128.nc
erf.rrtmgp_cloud_optics_sw = rrtmgp-cloud-optics-coeffs-sw.nc
erf.rrtmgp_cloud_optics_lw = rrtmgp-cloud-optics-coeffs-lw.nc
erf.rad_freq_in_steps      = 600 # number of time steps before running radiation model
erf.rad_t_sfc              = 300.0
