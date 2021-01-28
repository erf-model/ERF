Input Files and Controls
------------------------

The input file specified on the command line is a free-format text file, one entry per row, that specifies input data processed by the AMReX ``ParmParse`` module.

This file needs to be specified along with the executable as an ``argv`` option, for example:


::

	mpirun -np 64 ./ERF2d.xxx.yyy.ex inputs

Also, any entry that can be specified in the inputs file can also be specified on the command line; values specified on the command line override values in the inputs file, e.g.:

::

	mpirun -np 64 ./ERF2d.gnu.DEBUG.MPI.ex inputs amr.restart=sod_x_chk0030 erf.riemann_solver=3

The available options are divided into groups: those that control primarily AMReX are prefaced with ``amr.`` while those that are specific to ERF are prefaced with ``erf.``.

A typical input file looks something like the example below; 
a full list of ERF-specific input parameters are in ``Source/_cpp_parameters``. 
These parameters, once read, are available in the ``ERF`` object for use from C++ and are also copied to the module ``prob_params_module`` for use in FORTRAN. 

::

    # ------------------  INPUTS TO MAIN PROGRAM  -------------------
    #absolute stop time for the simulation
    stop_time = 6 

    #maximum number of time steps at base AMR level
    max_step = 30 
    # ---------------------------------------------------------------
    
    #------------------------
    # PROBLEM SIZE & GEOMETRY
    # -----------------------

    #flag for periodicity (here x direction is periodic)
    geometry.is_periodic = 1 0 0  
    
    #0 => cart, 1 => RZ  2=>spherical
    geometry.coord_sys   = 0      

    #coordinates of domain's lower corner
    geometry.prob_lo     =   -0.3     0.0   0.0     

    #coordinates of domain's upper corner
    geometry.prob_hi     =    0.3     0.3   0.15  

    #number of cells along each direction at base level (note: dx=dy=dz)
    amr.n_cell           =    128     64    32   
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    ERF specific inputs
    # ---------------------------------------------------------------

    # 0: Collela, Glaz and Ferguson (default)
    # 1: Collela and Glaz  
    # 2: HLLC
    erf.riemann_solver    = 0     

    # >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
    # Interior, UserBC, Symmetry, SlipWall, NoSlipWall
    # >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<

    #boundary condition at the lower face of each coordinate direction
    erf.lo_bc       =  "Interior"  "UserBC"  "SlipWall"        
    
    #boundary condition at the upper face of each coordinate direction
    erf.hi_bc       =  "Interior"  "UserBC"  "SlipWall"          
    
    #------------------------
    # TIME STEP CONTROL
    #------------------------

    erf.cfl            = 0.5     # cfl number for hyperbolic system
    erf.init_shrink    = 0.3     # first timestep is scaled by this factor
    erf.change_max     = 1.1     # maximum factor by which timestep can increase
    erf.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt

    #------------------------
    # WHICH PHYSICS
    #------------------------
    
    erf.do_hydro = 1               # enable hyperbolic term
    erf.do_mol_AD = 1              # use method of lines (MOL)
    erf.do_react = 0               # enable chemical reactions
    erf.ppm_type = 2               # piecewise parabolic reconstruction type
    erf.allow_negative_energy = 0  # flag to allow negative internal energy
    erf.diffuse_temp = 0           # enable thermal diffusion
    erf.diffuse_vel  = 0           # enable viscous diffusion
    erf.diffuse_spec = 0           # enable species diffusion
    
    #------------------------
    # DIAGNOSTICS & VERBOSITY
    #------------------------
    
    # coarse time steps between computing integral of 
    # conserved variables in the  domain
    # these values should stabilize at steady state
    erf.sum_interval = 1       

    erf.v            = 1        # verbosity in ERF cpp files
    amr.v              = 1        # verbosity in Amr.cpp
    #amr.grid_log       = grdlog  # name of grid logging file
    # ---------------------------------------------------------------
    
    # ---------------------------------------------------------------
    AMR specific inputs
    # ---------------------------------------------------------------
    
    #------------------------
    # REFINEMENT / REGRIDDING 
    #------------------------
    
    amr.max_level       = 2       # maximum level number allowed
    amr.ref_ratio       = 2 2 2 2 # refinement ratio across levels
    amr.regrid_int      = 2 2 2 2 # how often to regrid
    amr.blocking_factor = 8       # block factor in grid generation
    amr.max_grid_size   = 64      # maximum number of cells per box along x,y,z
    
    #specify species name as flame tracer for 
    #refinement purposes
    erf.flame_trac_name = HO2
    
    #------------------------
    # CHECKPOINT FILES
    #------------------------

    amr.checkpoint_files_output = 1
    amr.check_file              = chk    # root name of checkpoint/restart file
    amr.check_int               = 500    # number of timesteps between checkpoints
    
    #------------------------
    # PLOTFILES
    #------------------------
    
    amr.plot_files_output = 1
    amr.plot_file         = plt     # root name of plotfile
    amr.plot_int          = 100     # number of timesteps between plotfiles

    #pick which all derived variables to plot
    amr.derive_plot_vars  = pressure x_velocity y_velocity
    
    # probin filename that has tagging and other namelists
    amr.probin_file = probin 
    # ---------------------------------------------------------------
    
