# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Station time-series output (erf.station_names).
#
# The Straka density current, refined over the lower middle of the domain, with
# four stations chosen to exercise every path the sampler has:
#
#   Center  two locations inside the refined region, so the values must come
#           from level 1, and two heights, so the vertical interpolation runs
#   Edge    inside the outer half cell of the non-periodic x boundary, where
#           the horizontal stencil collapses to the edge cell
#   Wrap    inside the outer half cell of the periodic y boundary, where the
#           stencil reaches across the periodic image
#   Surface a 2D diagnostic, which has no height and is filled by the 2D
#           plotfile path rather than the 3D one
#
# The refined box stops at z = 3200 of a 6400 m domain, so level 1 covers only
# the lower half of its own vertical extent.  Both Center heights are inside
# that half, so level 1 supplies them: a level qualifies when it covers the
# column from the bottom of the domain up through the cells the vertical
# interpolation reads, not the whole column.  Run with erf.v = 1, the sampler
# prints the level it chose for each station, which is what makes that visible.
#
# Run by add_test_box_parity (the series must not depend on the decomposition)
# and by add_test_restart_parity (a restart must append a continuous series).

erf.prob_name = "Density Current"

erf.init_type = Isentropic

max_step   = 10
stop_time  = 900.0

erf.buoyancy_type = 1

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -12800.   0.    0.
geometry.prob_hi     =  12800. 100. 6400.
amr.n_cell           =  256      4    64     # dx=dy=dz=100 m, Straka et al 1993

geometry.is_periodic = 0 1 0

xlo.type = "Symmetry"
xhi.type = "Outflow"

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt       = 1.0
erf.fixed_fast_dt  = 0.25
erf.vert_implicit  = false

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = -1
erf.v              = 1
amr.v              = 1

# REFINEMENT / REGRIDDING
amr.max_level        = 1
amr.ref_ratio_vect   = 2 1 2
erf.regrid_int       = 2
erf.refinement_indicators = box1
erf.box1.max_level   = 1
erf.box1.in_box_lo   = -4000.   0.    0.
erf.box1.in_box_hi   =  4000. 100. 3200.

# CHECKPOINT FILES
erf.check_file      = chk
erf.check_int       = -1

# PLOTFILES
erf.plot_file_1     = plt
erf.plot_int_1      = 10
erf.plot_vars_1     = density x_velocity y_velocity z_velocity theta

# SOLVER CHOICE
erf.use_gravity = true
erf.use_coriolis = false
erf.les_type = "None"
erf.molec_diff_type = "ConstantAlpha"
erf.rho0_trans = 1.0
erf.dynamic_viscosity = 75.0
erf.alpha_T = 75.0

prob.T_pert = -15.0
prob.x_c = 0.0
prob.z_c = 3000.0
prob.x_r = 4000.0
prob.z_r = 2000.0
prob.T_0 = 300.0

# STATION TIME SERIES
erf.station_names = Center Edge Wrap Surface

erf.Center.field      = theta magvel
erf.Center.x          = -1000.0 1000.0
erf.Center.y          =     50.0   50.0
erf.Center.height_agl =    500.0 3000.0

erf.Edge.field        = theta
erf.Edge.x            = -12790.0
erf.Edge.y            =     50.0
erf.Edge.height_agl   =    500.0

erf.Wrap.field        = theta
erf.Wrap.x            =   -500.0
erf.Wrap.y            =      5.0
erf.Wrap.height_agl   =    500.0

erf.Surface.field     = z_surf
erf.Surface.x         =      0.0
erf.Surface.y         =     50.0

# Required: naming stations with no cadence stops the run, as it does for the
# line and plane samplers.  Every step here, so that a ten-step run gives the
# parity tests a series with something in it.
erf.station_sampling_interval = 1

# Flush more often than the 100-step default, so that a ten-step run exercises
# the flush path and, on a restart, the header check that fires with the first
# flush of the restarted run.
erf.station_buffer_steps = 2
