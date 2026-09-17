# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Box, rank and tiling parity of the turbulence closures and PBL schemes.
#
# The unstable, perturbed ABL of ABL_MRF_Tiling (32 x 32 x 32, surface layer, geostrophic
# wind), run twice by add_test_box_parity: once on one box, on one rank, without tiling,
# and once on four 16 x 16 x 32 boxes on two ranks with 8 x 8 tiles.  The plotfiles must
# agree.  The CTest entries select the closure (Deardorff, k-eqn, MYNN25, MYNNEDMF, MYJ,
# native SHOC), the moisture model, the mesh and the integrator on the command line.
# The boxes are never split in z, so the column solves (implicit diffusion, substep) apply.

erf.prob_name = "ABL"

max_step = 10

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 6400  6400  1024
amr.n_cell           =   32    32    32
# Four 16x16x32 boxes: full columns (MRF requires it), wider than a tile in y
amr.max_grid_size_x  = 16
amr.max_grid_size_y  = 16
amr.max_grid_size_z  = 32
geometry.is_periodic = 1 1 0

# SURFACE LAYER (heated, so the countergradient terms are active)
zlo.type = "surface_layer"
erf.most.z0             = 0.1
erf.most.surf_temp_flux = 0.24
# No PBL height by default; the entries whose scheme needs one (k-eqn length cap,
# MYNN25, MYJ) ask for MYNN25 on the command line
erf.most.pblh_calc      = "None"
# Local (per-column) averaging, so u_star, t_star and Olen vary with the
# perturbations instead of being one plane-averaged number in every column
erf.most.average_policy = 1

zhi.type       = "SlipWall"
zhi.theta_grad = 0.003

# INITIALIZATION
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "sounding_dry"

# Deterministic (position-only) perturbations so the columns differ and the
# initial state does not depend on the decomposition.  A field that is uniform
# across a plane would compare equal whatever the decomposition did to it, so
# the perturbations are what gives the parity check something to catch.  The
# sounding puts a 6 K inversion at 150-250 m and the wind is 5 m/s, so the
# bulk-Richardson crossing that sets the PBL height of the schemes that need
# one also lies inside the perturbed layer (z <= 200 m).
prob.pert_ref_height = 200.0
prob.pert_deltaU     = 1.0
prob.pert_deltaV     = 1.0
prob.pert_periods_U  = 3.0
prob.pert_periods_V  = 2.0
prob.pert_deltaT     = 0.5
prob.pert_periods_T  = 2.0

# TIME STEP CONTROL (dx = dy = 200 m, dz = 32 m: at dx = 100 m with ratio 4
# the acoustic substep is unstable and the run goes to NaN by step 5)
erf.fixed_dt           = 1.0
erf.fixed_mri_dt_ratio = 6

# DIAGNOSTICS & VERBOSITY
erf.sum_interval = 1
erf.v            = 1
amr.v            = 1

# REFINEMENT / REGRIDDING
amr.max_level = 0

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES (the two runs write the same prefix; RunBoxParity.cmake keeps them
# apart by running each in its own directory, one_box/ and split/)
erf.plot_file_1   = plt
erf.plot_int_1    = 10
erf.plot_vars_1   = density x_velocity y_velocity z_velocity pressure theta KE Kmv Khv

# SOLVER CHOICE
erf.molec_diff_type = "None"
erf.use_gravity     = true

erf.use_coriolis           = true
erf.latitude               = 45.0
erf.rotational_time_period = 86455.2516813368

erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind    = 5.0 0.0 0.0

# TURBULENCE CLOSURE
erf.les_type = "None"
erf.pbl_type = "None"
