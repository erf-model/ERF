# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 20

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 200.   40.  400.
amr.n_cell           =  40     8     79

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TERRRAIN GRID TYPE
erf.use_terrain = true      # enable terrain stencils
erf.terrain_type = Moving   # moving terrain
erf.terrain_smoothing = 2   # Sullivan 2004 approach

erf.no_substepping     = 1
erf.fixed_dt           = 0.005

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = -1        # timesteps between computing mass
erf.v              =  1        # verbosity in ERF.cpp
amr.v              =  1        # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk     # root name of checkpoint file
erf.check_int       = -1      # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt     # prefix of plotfile name
erf.plot_int_1      = 10      # number of timesteps between plotfiles
#erf.plot_vars_1     = x_velocity z_velocity pert_pres
    
erf.plot_vars_1     = density x_velocity y_velocity z_velocity pressure theta pres_hse dens_hse z_phys detJ 

# SOLVER CHOICE
erf.use_gravity = true
erf.use_coriolis = false
erf.use_rayleigh_damping = false

erf.buoyancy_type = 1

erf.dycore_horiz_adv_type  = Centered_2nd
erf.dycore_vert_adv_type   = Centered_2nd
erf.dryscal_horiz_adv_type = Centered_2nd
erf.dryscal_vert_adv_type  = Centered_2nd

# MULTILEVEL
amr.max_level = 0

#erf.terrain_z_levels = 0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.04150247509127 12.126229880712515 13.255976158466375 14.43260970297217 15.658076451850373 16.934403103949077 18.263700471134833 19.648166969191106 21.09009225359766 22.59186100620377 24.15595687905777 25.78496660191525 27.481584260219066 29.248615750626044 31.088983421449132 33.00573090568933 35.002028154650475 37.0811766804614 39.246615016175845 41.501924402479965 43.85083471041228 46.29723060989119 48.84515799425161 51.49883067141582 54.26263733276442 57.14114881123273 60.13912564063612 63.26152592872602 66.51351355699737 69.90046672080815 73.42798682393531 77.10190774227647 80.92830547201837 84.91350817822818 89.06410666048615 93.38696525286788 97.88923317630304 102.57835636208479 107.46208976608435 112.54851019403581 117.84602965910207 123.36340929381277 129.10977383938211 135.09462673636932 141.32786584163867 147.8197997976124 154.58116508088813 161.62314375841598 168.9573819806012 176.59600924191614 184.5516584408753 192.83748677254914 201.46719748816884 210.45506255780842 219.81594627362443 229.56532983268838 239.71933694006762 250.29476047449748 261.30909026074465 272.7805419945926 284.72808736828534 297.17148544625155 310.13131534299936 323.629010257224 337.6868929184159 352.3282125045899 367.5771830921924 383.4590237017744 400.00

erf.molec_diff_type = "None"
erf.les_type   = "None"

# PROBLEM PARAMETERS (optional)
prob.Ampl = 0.16
