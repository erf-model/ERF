# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 100000000
stop_time = 0.0018336339443081453
max_step = 100

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 1 1 1
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     =  -1.0 -1.0 -1.0
geometry.prob_hi     =   1.0  1.0  1.0
# use with single level
amr.n_cell           =  64    64    64
# use with 1 level of refinement
#amr.n_cell           =  128   128   128

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
erf.lo_bc       =  "Interior"  "Interior"  "Interior"
erf.hi_bc       =  "Interior"  "Interior"  "Interior"

# WHICH PHYSICS
erf.do_mol=1
erf.do_hydro = 1
erf.diffuse_vel = 1
erf.diffuse_temp = 1
erf.do_react = 0
erf.do_grav = 0
erf.use_explicit_filter=1
erf.les_filter_type=2
erf.les_filter_fgr=2

# TIME STEP CONTROL
erf.cfl            = 0.9     # cfl number for hyperbolic system
erf.init_shrink    = 0.3     # scale back initial timestep
erf.change_max     = 1.1     # max time step growth
erf.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in Castro.cpp
amr.v                = 1       # verbosity in Amr.cpp
amr.data_log         = datlog
#amr.grid_log        = grdlog  # name of grid logging file

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed
#amr.max_level       = 1       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2 2 2 2 # how often to regrid
amr.blocking_factor = 4       # block factor in grid generation
amr.max_grid_size   = 16
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est

# CHECKPOINT FILES
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 100        # number of timesteps between checkpoints

# PLOTFILES
amr.plot_file       = plt        # root name of plotfile
amr.plot_int        = 100        # number of timesteps between plotfiles
amr.plot_vars  =  density Temp
amr.derive_plot_vars = x_velocity y_velocity z_velocity magvel magvort pressure

# PROBLEM PARAMETERS
prob.reynolds = 1600.0
prob.mach = 0.1
prob.prandtl = 0.71

# EB
eb2.geom_type = "all_regular"
ebd.boundary_grad_stencil_type = 0
