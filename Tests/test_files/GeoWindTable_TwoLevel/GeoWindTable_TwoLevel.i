# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# A fluid at rest under a geostrophic wind read from a table (erf.abl_geo_wind_table),
# Ug = 2 + z/100 m/s, Vg = 0, on two levels with the fine level refined in z as well.
# After one step v = f dt Ug(z) in every cell of every level (f = coriolis factor *
# sin(latitude)); Tests/check_geowind_profile.py checks it.
erf.prob_name = "ABL"

max_step = 1
amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 800  800  1000
amr.n_cell           =   8    8    16
amr.max_grid_size    =   8    8    16
geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# INITIALIZATION: uniform state at rest
erf.init_type = "Uniform"
prob.rho_0 = 1.0
prob.T_0   = 300.0
prob.U_0   = 0.0
prob.V_0   = 0.0
prob.W_0   = 0.0

# TIME STEP CONTROL
erf.fixed_dt           = 1.0
erf.fixed_mri_dt_ratio = 4

# DIAGNOSTICS & VERBOSITY
erf.v = 0
amr.v = 0

# REFINEMENT: level 1 covers the whole horizontal extent (no coarse/fine side faces, which
# would add their own small velocity changes) and is refined in z as well
amr.max_level        = 1
amr.ref_ratio_vect   = 2 2 2
erf.refinement_indicators = box1
erf.box1.max_level   = 1
erf.box1.in_box_lo   =   0.    0.    0.
erf.box1.in_box_hi   = 800.  800. 1000.

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 1
erf.plot_vars_1 = density x_velocity y_velocity z_velocity

# SOLVER CHOICE
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.use_gravity     = false

erf.use_coriolis           = true
erf.latitude               = 90.0
erf.rotational_time_period = 125663.7061435917   # coriolis factor 1e-4 1/s

erf.abl_driver_type    = "GeostrophicWind"
erf.abl_geo_wind_table = "geo_wind_table"
