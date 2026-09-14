# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Flat periodic column with a stretched vertical grid (dz0 = 10 m, ratio 1.1,
# 40 levels) and the MOST surface layer with erf.most.zref unset.
#
# The default MOST query height (10 m) is exactly the top face of the first
# cell, and the uniform (prob_hi - prob_lo)/nz spacing is 110.6 m, so this deck
# exercises both reference-height defects fixed together:
#   - terrain-fitted lookups rejected a height on a face ("zref not found")
#   - the no-terrain default took half of the uniform spacing as zref
# run_most_zref.py runs it in several mesh configurations and checks u* against
# the log law and against a uniform 10 m column.
erf.prob_name = "ABL"
max_step = 10

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 400 400 4425.925556817605
amr.n_cell           = 4 4 40
amr.max_grid_size_z  = 512
geometry.is_periodic = 1 1 0

erf.initial_dz            = 10.0
erf.grid_stretching_ratio = 1.1

zlo.type = "surface_layer"
erf.most.z0 = 0.1
erf.most.surf_temp_flux = 0.0

zhi.type = "SlipWall"
zhi.theta_grad = 0.003

# TIME STEP CONTROL
# Compressible on purpose: the anelastic Poisson solve on a stretched or
# terrain-fitted mesh needs an FFT build, and the CI builds without FFT.
# The reference-height lookups do not depend on the dycore.
erf.anelastic          = 0
erf.fixed_dt           = 1.0
erf.fixed_mri_dt_ratio = 8

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.data_log     = hist.dat
erf.profile_int  = 1
erf.v            = 1
amr.v            = 1
amr.max_level    = 0

# CHECKPOINT & PLOTFILES
erf.check_int  = -1
erf.plot_int_1 = -1

# SOLVER CHOICE
erf.use_gravity     = true
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.pbl_type        = "MRF"

erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_most_zref"

erf.use_coriolis    = true
erf.latitude        = 45.0
erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind    = 15.0 0.0 0.0
