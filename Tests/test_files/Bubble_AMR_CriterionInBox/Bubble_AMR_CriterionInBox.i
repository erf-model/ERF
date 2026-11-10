# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# A field criterion restricted to a box, alongside a static box at the same level.
#   hot   tags theta > 300.5 but only for x >= 10 km, so only the +x half of the
#         bubble (centered at x = 10 km) is refined, and the box is not registered
#         as a static refined region
#   fixed is a static box; it must survive hot being read after it at every regrid
erf.prob_name = "Bubble"

max_step  = 10
stop_time = 3600.0

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 20000.0 20000.0 9600.0
amr.n_cell           = 32      32      16
geometry.is_periodic = 0 0 0
xlo.type = "SlipWall"
xhi.type = "SlipWall"
ylo.type = "SlipWall"
yhi.type = "SlipWall"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
erf.fixed_dt = 1.0

# DIAGNOSTICS & VERBOSITY
erf.v = 1
amr.v = 1

# REFINEMENT / REGRIDDING
amr.max_level       = 1
amr.ref_ratio       = 2
amr.n_error_buf     = 0
amr.grid_eff        = 1.0
amr.blocking_factor = 2
amr.max_grid_size   = 64
amr.refine_whole_domain_dir = 2
erf.regrid_int      = 5

erf.refinement_indicators = fixed hot

erf.fixed.max_level = 1
erf.fixed.in_box_lo = 0.    0.    6000.
erf.fixed.in_box_hi = 5000. 5000. 9600.

erf.hot.max_level     = 1
erf.hot.field_name    = theta
erf.hot.value_greater = 300.5
erf.hot.in_box_lo     = 10000.     0.    0.
erf.hot.in_box_hi     = 20000. 20000. 9600.

# CHECKPOINT FILES
erf.check_int = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 10
erf.plot_vars_1 = density rhotheta x_velocity y_velocity z_velocity theta pert_pres

# SOLVER CHOICES
erf.use_gravity = true

erf.dycore_horiz_adv_type  = "Upwind_3rd"
erf.dycore_vert_adv_type   = "Upwind_3rd"
erf.dryscal_horiz_adv_type = "Upwind_3rd"
erf.dryscal_vert_adv_type  = "Upwind_3rd"

# PHYSICS OPTIONS
erf.les_type        = "None"
erf.pbl_type        = "None"
erf.init_type       = Isentropic
erf.molec_diff_type = "None"

# PROBLEM PARAMETERS
prob.T_pert = 2.0
prob.x_c    = 10000.0
prob.y_c    = 10000.0
prob.z_c    =  2000.0
prob.x_r    =  2000.0
prob.y_r    =  2000.0
prob.z_r    =  2000.0
prob.T_0    =   300.0
prob.do_moist_bubble   = false
prob.T_pert_is_airtemp = false
