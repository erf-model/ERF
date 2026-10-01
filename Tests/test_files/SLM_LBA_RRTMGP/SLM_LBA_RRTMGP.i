# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 21600 # 6 hours

erf.prob_name = "SLM"

start_datetime = "1999-02-23 11:30:00"

amrex.fpe_trap_invalid = 0

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  16000. 16000. 25000.
amr.n_cell           =  32     32    250     # dx=dy=dz=500 m

geometry.is_periodic = 1 1 0

amrex.use_gpu_aware_mpi=1
amrex.the_arena_init_size = 1e9

erf.fix_random_seed = 1
erf.use_Ri_correction = 0

erf.mg_v = 1
erf.terrain_type = StaticFittedMesh
erf.flat_terrain = true
#erf.anelastic = 1
#erf.use_fft = 1

# MOST BOUNDARY (DEFAULT IS ADIABATIC FOR THETA)
zlo.type = "surface_layer"
erf.surface_layer.flux_type = "custom"
erf.most.ustar  = 0.  # ustar
erf.most.tstar  = 0. # theta flux 260W/m2
erf.most.qstar  = 0. # qv    flux 536 W/m2
erf.most.z0     = 0.1
erf.most.zref   = 50.0

#  -- Surface forcing for MOST
#erf.most.use_sfc_fluxes = 0 #1
#erf.most.sfc_file = "sfc_cass"

# NOTE: This should have a qv grad too (use hoextrapcc?!)
zhi.type = "SlipWall"
#zhi.theta_grad = 0.00365

# TIME STEP CONTROL
#erf.fixed_dt           = 5.0 # (anelastic dt) fixed time step depending on grid resolution
erf.fixed_dt            = 1.0 # fixed time step depending on grid resolution
erf.substepping_cfl = 0.5
#erf.fixed_mri_dt_ratio = 4

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp
erf.data_log       = "surf" "mean" "flux" "subgrid" "forcing"
#erf.profile_int    = 12     # (every minute with dt = 0.075)
erf.profile_int    = 60     # (every minute with dt = 0.075)

#erf.destag_profiles = false

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk   # root name of checkpoint file
erf.check_int       = -1    # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt   # prefix of plotfile name
erf.plot_int_1      = 300   # number of timesteps between plotfiles
#erf.plot_int_1      = 12   # number of timesteps between plotfiles
erf.plot_vars_1     = z_phys density rhoadv_0 rhotheta x_velocity y_velocity z_velocity pressure temp theta qt qp qv qc qi qsat scalar qsrc_lw qsrc_sw
#erf.plot_vars_1     = ttend qtend wsub tnudge qnudge unudge wnudge

erf.plotfile_type = "amrex"
erf.plot_lsm = true

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = true

#erf.use_coriolis    = false
erf.use_coriolis    = true
erf.coriolis_3d     = false
erf.latitude        = 14.982176712702886  # f = 0.376e-4 1/s

#erf.dycore_horiz_adv_type    = Upwind_3rd
#erf.dycore_vert_adv_type     = Upwind_3rd
#erf.dryscal_horiz_adv_type   = Upwind_3rd
#erf.dryscal_vert_adv_type    = Upwind_3rd
#erf.moistscal_horiz_adv_type = WENO5
#erf.moistscal_vert_adv_type  = WENO5

erf.dycore_horiz_adv_type    = Upwind_3rd
erf.dycore_vert_adv_type     = Upwind_3rd
erf.dryscal_horiz_adv_type   = Upwind_3rd_SL
erf.dryscal_vert_adv_type    = Upwind_3rd_SL
erf.moistscal_horiz_adv_type = Upwind_3rd_SL
erf.moistscal_vert_adv_type  = Upwind_3rd_SL


erf.moisture_model  = "Kessler_NoRain"
#erf.moisture_model  = "SatAdj"
#erf.buoyancy_type   = 4
erf.buoyancy_type   = 1

erf.molec_diff_type = "None"

erf.les_type        = "Smagorinsky"
erf.Cs              = 0.15

#erf.Pr_t      = 1.0
#erf.Sc_t      = 1.0

erf.Pr_t      = 0.33333333333333
erf.Sc_t      = 0.33333333333333

erf.init_type = "input_sounding"
erf.sounding_type = Ideal

erf.input_sounding_file = "snd_lba"
erf.input_sounding_time = 0.0

# nudging
#erf.nudging_from_input_sounding = true
erf.nudging_from_input_sounding = false

erf.tau_nudging = 3600.0 # tautqls

# set below 4 to 30K to bypass t and q nudging
erf.nudging_t_z1 = 30000.0
erf.nudging_t_z2 = 30000.0
erf.nudging_q_z1 = 30000.0
erf.nudging_q_z2 = 30000.0

# large scale forcing
erf.large_scale_forcing = false

erf.land_surface_model = "SLM"

# ------------------------
# SLM configuration
# ------------------------
# number of soil layers
slm.nsoil = 9
# soil thickness (in m) for each layer (from top layer to bottom)
slm.soil_dz = 0.01 0.02 0.04 0.08 0.16 0.32 0.37 0.50 1.00

# initial values used over entire grid
slm.landtype0 = 2     # evergreen forest
slm.LAI0 = 6.        # default leaf area index
slm.clay0 = 13.0      # uniform clay content (%)
slm.sand0 = 17.0      # uniform sand content (%)
slm.sw0 = 0.6000000 0.6030075 0.6090226 0.6210526 0.6451128 0.6932331 0.7624060 0.8496241 1.0000000        # uniform soil wetness fraction
slm.st0 = 300.1500000 300.1374948 300.1124843 300.0624635 299.9624217 299.7623382 299.4747182 299.1120669 298.4868060             # uniform soil temperature (K)
slm.soiltnudging = true
slm.soilwnudging = true
slm.relax_hgt = 0 0 0 0 0 0 0 0 1
slm.tabs_s = 0.0      # surface temperature (K)
slm.t00 = 300.0       # constant offset for sstxy

slm.z0_soil = 0.0387  # baresoil roughness length
slm.mws_mx0 = 50.0    # maximum puddle water storage (mm)
slm.Rc_max  = 5000.0  # maximum stomatal resistance
slm.T_opt   = 298.0   # optimum temperature for transpiration
#slm.zref    = 25.0    # height of reference level (m)


#slm.rad_input_file = "LBA_240x240x250_default_rad_coszrs_combined.nc"

## Radiation
# RADIATION INPUTS
erf.radiation_model = "RRTMGP"
erf.rad.datalog = "radiation"
erf.profile_rad_int = 300
#erf.profile_rad_int = 120
erf.rad_t_sfc = 300
#erf.rad_freq_in_steps      = 120
erf.rad_freq_in_steps      = 300
erf.rad_do_subcol_sampling = false
erf.rad_write_fluxes       = false
erf.rad_orbital_year       = 1999
erf.rad_cons_lat           = -10.
erf.rad_cons_lon           = -60.
erf.co2vmr                 = 0.00036
erf.o3vmr                  = 2.9114E-08  2.9416E-08  2.9717E-08  3.0019E-08  3.0321E-08  3.0623E-08  3.0924E-08  3.1226E-08  3.1528E-08  3.1829E-08  3.2071E-08  3.2252E-08  3.2433E-08  3.2614E-08  3.2795E-08  3.2976E-08  3.3157E-08  3.3338E-08  3.3519E-08  3.3700E-08  3.3851E-08  3.3971E-08  3.4092E-08  3.4213E-08  3.4333E-08  3.4454E-08  3.4575E-08  3.4696E-08  3.4816E-08  3.4937E-08  3.5027E-08  3.5088E-08  3.5148E-08  3.5208E-08  3.5269E-08  3.5329E-08  3.5389E-08  3.5450E-08  3.5510E-08  3.5570E-08  3.5721E-08  3.5963E-08  3.6204E-08  3.6445E-08  3.6687E-08  3.6928E-08  3.7169E-08  3.7411E-08  3.7652E-08  3.7894E-08  3.8105E-08  3.8286E-08  3.8467E-08  3.8648E-08  3.8829E-08  3.9010E-08  3.9191E-08  3.9372E-08  3.9553E-08  3.9734E-08  3.9945E-08  4.0186E-08  4.0428E-08  4.0669E-08  4.0911E-08  4.1152E-08  4.1393E-08  4.1635E-08  4.1876E-08  4.2117E-08  4.2359E-08  4.2600E-08  4.2841E-08  4.3083E-08  4.3324E-08  4.3565E-08  4.3807E-08  4.4048E-08  4.4290E-08  4.4531E-08  4.4923E-08  4.5466E-08  4.6009E-08  4.6552E-08  4.7095E-08  4.7638E-08  4.8181E-08  4.8725E-08  4.9268E-08  4.9811E-08  5.0384E-08  5.0987E-08  5.1591E-08  5.2194E-08  5.2798E-08  5.3401E-08  5.4004E-08  5.4608E-08  5.5211E-08  5.5815E-08  5.6629E-08  5.7655E-08  5.8681E-08  5.9706E-08  6.0702E-08  6.1788E-08  6.2814E-08  6.3840E-08  6.4866E-08  6.5891E-08  6.6977E-08  6.8184E-08  6.9391E-08  7.0598E-08  7.1805E-08  7.3011E-08  7.4218E-08  7.5425E-08  7.6632E-08  7.7839E-08  7.9045E-08  8.0252E-08  8.1459E-08  8.2666E-08  8.3873E-08  8.5079E-08  8.6286E-08  8.7493E-08  8.8700E-08  8.9907E-08  9.1113E-08  9.2320E-08  9.3527E-08  9.4734E-08  9.5941E-08  9.7147E-08  9.8354E-08  9.9561E-08  1.0077E-07  1.0197E-07  1.0378E-07  1.0620E-07  1.0861E-07  1.1103E-07  1.1344E-07  1.1585E-07  1.1827E-07  1.2068E-07  1.2309E-07  1.2551E-07  1.2762E-07  1.2943E-07  1.3124E-07  1.3305E-07  1.3486E-07  1.3667E-07  1.3848E-07  1.4029E-07  1.4210E-07  1.4391E-07  1.4994E-07  1.6020E-07  1.7046E-07  1.8072E-07  1.9098E-07  2.0123E-07  2.1149E-07  2.2175E-07  2.3201E-07  2.4227E-07  2.5463E-07  2.6912E-07  2.8360E-07  2.9808E-07  3.1256E-07  3.2704E-07  3.4152E-07  3.5601E-07  3.7049E-07  3.8497E-07  4.0880E-07  4.4199E-07  4.7518E-07  5.0836E-07  5.4155E-07  5.7480E-07  6.0823E-07  6.4141E-07  6.7460E-07  7.0779E-07  7.4822E-07  7.9649E-07  8.4476E-07  8.9303E-07  9.4130E-07  9.8958E-07  1.0378E-06  1.0861E-06  1.1344E-06  1.1827E-06  1.2370E-06  1.2973E-06  1.3577E-06  1.4180E-06  1.4783E-06  1.5387E-06  1.5990E-06  1.6594E-06  1.7197E-06  1.7800E-06  1.8464E-06  1.9188E-06  1.9912E-06  2.0636E-06  2.1360E-06  2.2084E-06  2.2809E-06  2.3533E-06  2.4257E-06  2.4981E-06  2.5795E-06  2.6700E-06  2.7606E-06  2.8511E-06  2.9416E-06  3.0321E-06  3.1226E-06  3.2131E-06  3.3036E-06  3.3941E-06  3.4816E-06  3.5661E-06  3.6506E-06  3.7350E-06  3.8195E-06  3.9040E-06  3.9885E-06  4.0730E-06  4.1574E-06  4.2419E-06  4.3234E-06  4.4018E-06  4.4802E-06  4.5587E-06  4.6371E-06  4.7156E-06  4.7940E-06  4.8725E-06  4.9509E-06  5.0293E-06
erf.n2ovmr                 = 2.106e-7
erf.covmr                  = 1.5e-7
erf.ch4vmr                 = 3.068e-6
erf.o2vmr                  = 0.209
erf.n2vmr                  = 0.7906
erf.nswbands               = 14
erf.nlwbands               = 16
erf.nswgpts                = 112 #224
erf.nlwgpts                = 128 #256
erf.rrtmgp_file_path       = ./
#erf.rrtmgp_coeffs_sw       = rrtmgp-data-sw-g224-2018-12-04.nc
#erf.rrtmgp_coeffs_lw       = rrtmgp-data-lw-g256-2018-12-04.nc
erf.rrtmgp_coeffs_sw       = rrtmgp-gas-sw-g112.nc
erf.rrtmgp_coeffs_lw       = rrtmgp-gas-lw-g128.nc
erf.rrtmgp_cloud_optics_sw = rrtmgp-cloud-optics-coeffs-sw.nc
erf.rrtmgp_cloud_optics_lw = rrtmgp-cloud-optics-coeffs-lw.nc
erf.rad_day0               = 0 #920023400 #919596600
# Higher values of perturbations lead to instability
# Instability seems to be coming from BC
prob.U_0_Pert_Mag = 0.00
prob.V_0_Pert_Mag = 0.00
prob.W_0_Pert_Mag = 0.0

prob.pert_ref_height = 1000.0
prob.T_0_Pert_Mag    = 0.1
prob.qv_0_Pert_Mag   = 0.0

prob.custom_TKE      = false
#prob.custom_TKE      = true
