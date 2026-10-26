# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Anelastic flow in a closed box with a thin immersed wall (zero x-flux faces), so the
# projection goes through project_velocity_tb, with two corrector passes (erf.ncorr = 2).
# Tests/check_thin_body_gradp.py compares the stored pressure gradient (dpdx/dpdy/dpdz,
# the cell average of gradp) with centred differences of the stored pressure increment
# (pert_pres, pp_inc): both are built from the same solves, so they must agree.
erf.prob_name = "Flow in a Box"

max_step = 2
amrex.fpe_trap_invalid = 1

erf.anelastic = 1
erf.ncorr     = 2
erf.init_type = Uniform

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo = -1. -1.  0.
geometry.prob_hi =  1.  1.  1.
amr.n_cell       = 16  16   8
amr.max_grid_size = 8 8 8

geometry.is_periodic = 0 0 0

xlo.type = "SlipWall"
xhi.type = "SlipWall"
ylo.type = "SlipWall"
yhi.type = "SlipWall"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

xlo.theta = 285.
xhi.theta = 285.
ylo.theta = 285.
yhi.theta = 285.
zlo.theta = 299.
zhi.theta = 280.

# Thin wall: x-faces i = 8, j = 6..9, k = 0..3
erf.zero_xflux_faces = (8,6,0) (8,7,0) (8,8,0) (8,9,0) (8,6,1) (8,7,1) (8,8,1) (8,9,1) (8,6,2) (8,7,2) (8,8,2) (8,9,2) (8,6,3) (8,7,3) (8,8,3) (8,9,3)

# TIME STEP CONTROL
erf.fixed_dt = 0.025

# DIAGNOSTICS & VERBOSITY
erf.v  = 0
amr.v  = 0

# REFINEMENT / REGRIDDING
amr.max_level = 0

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 2
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta pert_pres dpdx dpdy dpdz

# SOLVER CHOICE
erf.use_gravity = true
erf.molec_diff_type = "None"
erf.vert_implicit   = false
erf.les_type = "Smagorinsky"
erf.Cs       = 0.1
erf.Pr_t     = 0.33333333333333

# PROBLEM PARAMETERS
prob.rho_0        = 1.0
prob.T_0          = 285.
prob.T_0_Pert_Mag = 0.1
