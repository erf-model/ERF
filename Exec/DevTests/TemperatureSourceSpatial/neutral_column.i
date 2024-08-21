# ------------------  INPUTS TO MAIN PROGRAM  -------------------
stop_time = 999.9
max_step = 20

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent    =   640.  640.  160.
amr.n_cell              =   64    64    16

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# INITIALIZATION
erf.init_type = input_sounding
erf.init_sounding_ideal = true

# TIME STEP CONTROL
erf.fixed_dt                    = 0.1

# DIAGNOSTICS & VERBOSITY
amr.v               = 1     # verbosity in Amr.cpp
erf.v               = 1     # verbosity in ERF.cpp -- needs to be 1 to write out data_log files
erf.sum_interval    = 1     # timesteps between computing mass

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = -1         # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt       # prefix of plotfile name
erf.plot_int_1      = 5         # number of timesteps between plotfiles (DEBUG)
erf.plot_vars_1     = density rhotheta x_velocity y_velocity z_velocity pressure theta temp rhoQ1 qt qv qc qi qp qsat qrain rain_accum snow_accum graup_accum


# SOLVER CHOICES
erf.use_gravity = true
erf.use_coriolis = false

erf.moisture_model  = "SAM"
erf.buoyancy_type   = 1

erf.molec_diff_type = "None"
erf.les_type        = "Smagorinsky"
erf.Cs              = 0.17

erf.add_custom_rhotheta_forcing = true
erf.add_custom_moisture_forcing = true
erf.spatial_rhotheta_forcing = true
erf.spatial_moisture_forcing = true
erf.custom_forcing_uses_primitive_vars = true
prob.advection_heating_rate   = 2.3148E-4
prob.advection_moisture_rate  = 1.2E-8

