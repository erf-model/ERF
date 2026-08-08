# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Minimal TwoStream radiation test with prognostic cloud fraction
# Tests Phase 14A bugfix: TwoStream MultiFab vector allocation
# No LSM present; uses fallback scalar surface properties

erf.prob_name = "TwoStream Prognostic Cloud Fraction"

max_step = 5

amrex.fpe_trap_invalid = 1

erf.anelastic = 0

erf.init_type = Uniform

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo = -1. -1.  0.
geometry.prob_hi =  1.  1.  1.

#coarse
amr.n_cell       = 32  32  16

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
zlo.theta = 290.
zhi.theta = 280.

# TIME STEP CONTROL
erf.fixed_dt      = 0.01

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1
erf.v              = 1
amr.v              = 1

# REFINEMENT / REGRIDDING
amr.max_level       = 0

# CHECKPOINT FILES
erf.check_file      = chk
erf.check_int       = 1000

# PLOTFILES
erf.plot_file_1     = plt
erf.plot_int_1      = 100
erf.plot_vars_1     = density x_velocity y_velocity z_velocity pressure temp theta

erf.use_gravity = true

# SOLVER CHOICE
erf.molec_diff_type  = "None"
erf.alpha_T = 0.0
erf.alpha_C = 0.0

erf.les_type  = "Smagorinsky"
erf.Cs        = 0.1
erf.Pr_t      = 0.33333333333333

# RADIATION
erf.radiation = true
erf.radiation.rad_type = TwoStream
erf.radiation.surface_albedo_sw = 0.3
erf.radiation.surface_emissivity_lw = 0.99
erf.radiation.surface_temp_k = 290.0
erf.radiation.cloud_fraction_prog_enable = true

# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.T_0          = 285.
prob.T_0_Pert_Mag = 0.1
