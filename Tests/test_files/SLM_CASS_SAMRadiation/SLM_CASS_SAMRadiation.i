# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#stop_time = 3600    #  1 hour
#stop_time = 4800    #  1.5 hours
#stop_time = 7200    #  2 hours
#stop_time = 43200   #  12 hours
#max_step = 43200
max_step = 34500

erf.prob_name = "SLM"

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  1600. 1600. 16000.
amr.n_cell           =  32    32    320    # dx=dy=dz=50 m

geometry.is_periodic = 1 1 0

erf.anelastic = 1

erf.mg_v = 1
erf.terrain_type = StaticFittedMesh
erf.flat_terrain = true
erf.use_fft = true

amrex.use_gpu_aware_mpi=1

#erf.use_smag_stratification = false
#erf.mix_isotropic = false

# MOST BOUNDARY (DEFAULT IS ADIABATIC FOR THETA)
zlo.type = "surface_layer"
erf.surface_layer.flux_type = "custom"
erf.most.ustar  = 0.0
erf.most.tstar  = 0.0
erf.most.qstar  = 0.0
erf.most.z0     = 0.1
erf.most.zref   = 25.0

# NOTE: This should have a qv grad too (use hoextrapcc?!)
zhi.type = "SlipWall"
#zhi.theta_grad = 0.00365

# TIME STEP CONTROL
erf.fixed_dt           = 1.0 # fixed time step depending on grid resolution
erf.fixed_mri_dt_ratio = 4

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 0       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp
erf.data_log       = "surf" "mean" "flux" "subgrid" "forcing"
erf.profile_int    = 60     # (every minute with dt = 0.075)

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk     # root name of checkpoint file
erf.check_int       = -1    # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt     # prefix of plotfile name
erf.plot_int_1      = 300  #2300  # number of timesteps between plotfiles
erf.plot_vars_1     = density x_velocity y_velocity z_velocity pressure theta qt qp qv qc qsrc z_phys

erf.plotfile_type = "amrex"
erf.plot_lsm = true

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = true

erf.use_coriolis    = false
erf.coriolis_3d     = false
erf.latitude        = 14.982176712702886  # f = 0.376e-4 1/s

erf.dycore_horiz_adv_type    = Upwind_3rd
erf.dycore_vert_adv_type     = Upwind_3rd
erf.dryscal_horiz_adv_type   = Upwind_3rd
erf.dryscal_vert_adv_type    = Upwind_3rd
erf.moistscal_horiz_adv_type = WENO5
erf.moistscal_vert_adv_type  = WENO5


erf.moisture_model  = "Kessler_NoRain"
erf.buoyancy_type   = 4

erf.molec_diff_type = "None"

erf.les_type        = "Smagorinsky"
erf.Cs              = 0.15

#erf.Pr_t      = 0.33333333333333
#erf.Sc_t      = 0.33333333333333
erf.Pr_t      = 1.0
erf.Sc_t      = 1.0

erf.init_type = "input_sounding"
erf.init_sounding_ideal = true
erf.sounding_type = Ideal

#erf.input_sounding_file = "sounding_cass"
erf.input_sounding_file = "sounding_cass_interpolated"
erf.input_sounding_time = 0.0

erf.nudging_from_input_sounding = true
erf.tau_nudging = 3600.0

erf.nudging_t_z1 = 5500.0
erf.nudging_t_z2 = 15900.0
erf.nudging_q_z1 = 5500.0
erf.nudging_q_z2 = 15900.0

erf.large_scale_forcing = true
erf.forcing_timescale = 3600.0 # tauls
erf.large_scale_forcing_file = "lsf_cass"

erf.land_surface_model = "SLM"

# ------------------------
# SLM configuration
# ------------------------
# number of soil layers
slm.nsoil = 9
# soil thickness (in m) for each layer (from top layer to bottom)
slm.soil_dz = 0.01 0.02 0.04 0.08 0.16 0.32 0.37 0.50 1.00

# initial values used over entire grid
slm.landtype0 = 10     # evergreen forest
slm.LAI0 = 2.0        # default leaf area index
slm.clay0 = 13.0      # uniform clay content (%)
slm.sand0 = 17.0      # uniform sand content (%)
slm.sw0 = 0.6000000 0.6030075 0.6090226 0.6210526 0.6451128 0.6932331 0.7624060 0.8496241 1.0000000        # uniform soil wetness fraction
slm.st0 = 300.1500000 300.1374948 300.1124843 300.0624635 299.9624217 299.7623382 299.4747182 299.1120669 298.4868060             # uniform soil temperature (K)
slm.relax_hgt = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 # soil relaxation function for nudging
slm.soiltnudging = true  # soil temperature nudging
slm.soilwnudging = true  # soil moisture nudging
slm.tausoil = 86400.0 # soil nudging timescale
slm.tabs_s = 300.0      # surface temperature (K)
slm.t00 = 300.0       # constant offset for sstxy

slm.z0_soil = 0.0387  # baresoil roughness length
slm.mws_mx0 = 50.0    # maximum puddle water storage (mm)
slm.Rc_max  = 5000.0  # maximum stomatal resistance
slm.T_opt   = 298.0   # optimum temperature for transpiration
slm.zref    = 25.0    # height of reference level (m)


slm.rad_input_file = "CASS_32x32x156_50m_50m_1s_rad_coszrs_combined.nc"


# Higher values of perturbations lead to instability
# Instability seems to be coming from BC
prob.U_0_Pert_Mag = 0.00
prob.V_0_Pert_Mag = 0.00
prob.W_0_Pert_Mag = 0.0

prob.pert_ref_height = 200.0
prob.T_0_Pert_Mag    = 0.1
prob.qv_0_Pert_Mag   = 0.0

prob.custom_TKE      = false
#prob.custom_TKE      = true
