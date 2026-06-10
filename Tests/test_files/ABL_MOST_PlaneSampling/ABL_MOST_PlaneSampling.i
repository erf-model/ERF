# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Regression test for plane_sampling_vars in PlaneSampler.
# Based on ABL_MOST; adds a horizontal plane sample at mid-domain
# height with a non-default variable list to exercise the new feature.
erf.prob_name = "ABL"

max_step = 10

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  1024     1024    1024
amr.n_cell           =    64       64      64

geometry.is_periodic = 1 1 0

# MOST BOUNDARY (DEFAULT IS ADIABATIC FOR THETA)
zlo.type = "surface_layer"
erf.most.z0   = 0.1
erf.most.zref = 8.0

zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 0.1  # fixed time step depending on grid resolution

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = -1         # no checkpoint output

# PLOTFILES
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 10         # number of timesteps between plotfiles
erf.plot_vars_1     = density rhoadv_0 x_velocity y_velocity z_velocity pressure temp theta

# PLANE SAMPLING
erf.do_plane_sampling        = true
erf.plane_sampling_interval  = 5
erf.sample_plane_lo          =  48.0  48.0  512.0
erf.sample_plane_hi          = 976.0 976.0  512.0
erf.sample_plane_dir         = 2
erf.plane_sampling_vars      = theta x_velocity y_velocity magvel

# SOLVER CHOICE
erf.use_gravity = false

erf.molec_diff_type = "None"
erf.les_type  = "Deardorff"
erf.Ck        = 0.1
erf.sigma_k   = 1.0
erf.Ce        = 0.1
erf.theta_ref = 300.0

erf.init_type = Uniform

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.A_0 = 1.0

prob.U_0 = 10.0
prob.V_0 = 0.0
prob.W_0 = 0.0
prob.T_0 = 300.0

prob.U_0_Pert_Mag = 0.0
prob.V_0_Pert_Mag = 0.0
prob.W_0_Pert_Mag = 0.0
