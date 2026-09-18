# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Box parity of the TKE buoyancy source across z box boundaries: the convective
# Canonical RANS case (Canonical_RANS/Convective_ABL_Flat), compressible with
# explicit vertical diffusion, run on one box and on boxes split in z. The k
# source averages the theta-diffusion fluxes at the two z-faces of each cell, so
# the top cell of every box reads a face written for that box; the two runs must
# agree. amr.max_grid_size_z = 24 gives five boxes of 20 cells (BoxList::maxSize
# divides out the common factors of 2 first), so the splits are at k = 20, 40, 60
# and 80, that is z = 400, 800, 1200 and 1600 m: two below the 937 m inversion and
# two above it, with a non-zero heat flux at all four from the first step.
#
# Difference from the canonical deck: no acoustic substepping.  The implicit w
# solve of the substep closes every box's column with Dirichlet rows
# (ERF_MakeFastCoeffs.cpp), so it is only correct on unsplit columns; dt is then
# the acoustic limit.  The PBL height scans whole columns on split grids, so the
# k-eqn length cap it feeds (erf.rans_lscale_from_pblh) stays on.
erf.prob_name = "ABL"

erf.anelastic         = 0
erf.use_fft           = false
erf.substepping_type  = None
erf.vert_implicit     = false

max_step = 40

amrex.fpe_trap_invalid = 0

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_extent = 2560   2560   2000
amr.n_cell           =    8      8    100     # dz = 20 m, first cell centre 10 m
amr.max_grid_size_x  = 4                      # split run: 2 x 1 x 5 boxes
amr.max_grid_size_y  = 8
amr.max_grid_size_z  = 24
amr.blocking_factor  = 4

geometry.is_periodic = 1 1 0

# MOST BOUNDARY: prescribed surface heat flux
zlo.type                = "surface_layer"
erf.most.z0             = 0.16
erf.most.surf_temp_flux = 0.24   # [K m/s]
erf.most.pblh_calc      = "MYNN25"

zhi.type       = "SlipWall"
zhi.theta_grad = 0.003   # [K/m], matches the sounding above the inversion

# INITIALIZATION
erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"
erf.sounding_type       = "Ideal"

# TIME STEP CONTROL: acoustic limit dz / c = 0.058 s
erf.fixed_dt = 0.02

# DIAGNOSTICS & VERBOSITY
erf.v = 1
amr.v = 0

# REFINEMENT / REGRIDGING
amr.max_level = 0

# CHECKPOINT FILES
erf.check_int  = -1

# PLOTFILES
erf.plot_file_1 = plt
erf.plot_int_1  = 40
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta KE Kmv Khv

# SOLVER CHOICE
erf.dycore_horiz_adv_type  = "Upwind_5th"
erf.dycore_vert_adv_type   = "Upwind_3rd"
erf.dryscal_horiz_adv_type = "Upwind_5th"
erf.dryscal_vert_adv_type  = "Upwind_3rd"

erf.use_gravity  = true
erf.use_coriolis = true
erf.coriolis_3d  = false
erf.latitude     = 90.0
erf.rotational_time_period = 125663.7061435917   # f = 1e-4 1/s

erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind    = 10.0 0.0 0.0

# gravity-wave damping in the top 500 m
erf.rayleigh_damp_W     = true
erf.rayleigh_dampcoef   = 0.2
erf.rayleigh_zdamp      = 500.

erf.molec_diff_type = "None"

# TURBULENCE MODELING: one-equation k RANS
erf.les_type                      = "None"
erf.rans_type                     = "kEqn"
erf.dirichlet_k                   = true
erf.init_tke_from_ustar           = true
erf.rans_consistent_diffusivities = true
erf.rans_lscale_from_pblh         = true   # cap l_g at kappa * 0.1 * zi
erf.rans_lscale_min               = 1.0
erf.max_geom_lscale               = 100.0
erf.theta_ref                     = 300.0

# PROBLEM PARAMETERS: no perturbations
prob.pert_deltaU  = 0.0
prob.pert_deltaV  = 0.0
prob.T_0_Pert_Mag = 0.0
prob.U_0_Pert_Mag = 0.0
prob.V_0_Pert_Mag = 0.0
prob.W_0_Pert_Mag = 0.0
