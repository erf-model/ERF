# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Advection of a z-independent Gaussian blob by a uniform horizontal velocity
# on a *vertically stretched* mesh.
#
# The exact solution is pure horizontal translation, so the scalar must remain
# independent of z for all time.  The time-averaged mass flux that the acoustic
# substepping hands to scalar advection (avg_xmom/avg_ymom) has to carry the
# h_zeta = dz_k/dz metric weighting for that to hold -- without it each k-level
# advects at u/h_zeta(k) and the blob shears apart in the vertical.
#
erf.prob_name = "Scalar Advection/Diffusion"

max_step = 40

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent =  1     1     1
amr.n_cell           = 32    32    16

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# VERTICALLY STRETCHED MESH (geometric, ratio 1.15; dz varies by a factor of ~8
# so h_zeta = dz_k/dz spans [0.287, 2.337])
erf.terrain_type = StaticFittedMesh
erf.terrain_z_levels = 0.0000000000 0.0179476914 0.0385875364 0.0623233583 0.0896195534 \
                       0.1210101778 0.1571093958 0.1986234965 0.2463647124 0.3012671106 \
                       0.3644048685 0.4370132902 0.5205129751 0.6165376127 0.7269659460 \
                       0.8539585292 1.0000000000

# TIME STEP CONTROL
# Fixed dt (and a fixed number of acoustic substeps) so the test is bitwise
# reproducible; the acoustic substepping is what we are exercising here.
erf.substepping_type   = Implicit
erf.fixed_dt           = 0.0005
erf.fixed_mri_dt_ratio = 6

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = -1      # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = -1         # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 40         # number of timesteps between plotfiles
erf.plot_vars_1     = density rhoadv_0 scalar x_velocity y_velocity z_velocity pressure temp theta

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 0.0
erf.use_gravity = false

erf.les_type         = "None"
erf.molec_diff_type  = "None"

erf.init_type = "uniform"

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.T_0   = 1.0
prob.A_0   = 1.0
prob.U_0   = 10.0
prob.V_0   = 5.0
prob.rad_0 = 0.125
prob.uRef  = 0.0

# z-independent Gaussian passive scalar
prob.prob_type = 15
