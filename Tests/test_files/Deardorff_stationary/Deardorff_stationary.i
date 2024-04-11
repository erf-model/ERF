# ------------------  INPUTS TO MAIN PROGRAM  -------------------
stop_time = 999.9
max_step = 10

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent    =   125.    125.   1000.
amr.n_cell              =    16      16     128

geometry.is_periodic = 1 1 0

#zhi.type = "SlipWall"
#zhi.theta_grad = 0.0 # true neutral boundary layer
zhi.type = "NoSlipWall"
zhi.density = 1.0
zhi.theta = 290.0
zhi.velocity = 15 0 0 # to match input_sounding

#zlo.type = "SlipWall"
zlo.type = "NoSlipWall"
zlo.density = 1.0
zlo.theta = 290.0
zlo.velocity = 5 0 0 # to match input_sounding

# TIME STEP CONTROL
erf.fixed_dt                    = 0.05

# DIAGNOSTICS & VERBOSITY
amr.v               = 1     # verbosity in Amr.cpp
erf.v               = 1     # verbosity in ERF.cpp -- needs to be 1 to write out data_log files
erf.sum_interval    = 1     # timesteps between computing mass
erf.data_log        = scalars.hist h_avg_profiles1.hist h_avg_profiles2.hist h_avg_profiles3.hist
erf.profile_int     = 1

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = -1         # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt       # prefix of plotfile name
erf.plot_int_1      = 10        # number of timesteps between plotfiles (DEBUG)
erf.plot_vars_1     = density x_velocity y_velocity z_velocity pressure theta rhoKE #pres_hse dens_hse

# SOLVER CHOICES
erf.use_gravity = false
erf.use_coriolis = false

erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind = 0. 0. 0.  # no background pressure gradient

erf.molec_diff_type = "None"
erf.les_type = "Deardorff"
erf.Ck       = 0.1
erf.Ce       = 0.93
erf.Pr_t     = 0.3333
erf.theta_ref = 290.0 # used in buoyancy term
erf.KE_0  = 0.000656292002688172 # exact soln in uniform density field, e = Ck/Ce*(dUdz*delta)**2

# INITIAL PROFILES
erf.init_type = "input_sounding"
erf.input_sounding_file = "input_sounding" # with linear wind profile
