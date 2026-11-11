# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# A small stable inflow (from Exec/RegTests/TurbulentInflow/stable_5m_inflow_inputs_CPM) with
# inflow perturbations at the south face. The perturbation boxes (3 cells on a side) straddle the
# grid boxes of the split run, and the run is long enough for the boxes to refresh, so it checks
# that the perturbations do not depend on the domain decomposition (#4317).
# erf.perturbation_type is set by the test.
max_step = 30

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo =   0    0    0
geometry.prob_hi = 160  320  160
amr.n_cell       =  32   64   32

geometry.is_periodic = 1 0 0

ylo.type = "Inflow"
yhi.type = "Outflow"
ylo.dirichlet_file = "stable_inflow.txt"

zlo.type                    = "surface_layer"
erf.most.z0                 = 0.29
erf.most.k_arr_in           = 0 0
erf.most.pot_surf_temp      = 306.0
erf.most.surf_heating_rate  = -1.0

zhi.type = "SlipWall"
zhi.theta_grad = 0.0016666667

erf.sponge_type = "input_sponge"
erf.input_sponge_file = "stable_sponge.txt"
erf.sponge_strength = 2.5
erf.use_yhi_sponge_damping = true
erf.yhi_sponge_start = 280

# TIME STEP CONTROL
erf.cfl = 10.0
erf.substepping_cfl = 0.5

# DIAGNOSTICS & VERBOSITY
erf.v              = 0
amr.v              = 0

# REFINEMENT / REGRIDDING
amr.max_level       = 0

# CHECKPOINT FILES
erf.check_int       = -1

# PLOTFILES
erf.plot_file_1     = plt
erf.plot_int_1      = 30
erf.plot_vars_1     = density x_velocity y_velocity z_velocity theta

# SOLVER CHOICE
erf.dycore_horiz_adv_type  = "Upwind_3rd"
erf.dycore_vert_adv_type   = "Upwind_3rd"
erf.dryscal_horiz_adv_type = WENO3
erf.dryscal_vert_adv_type  = WENO3

erf.use_gravity = true
erf.use_coriolis = false

erf.molec_diff_type = "None"
erf.les_type        = "Deardorff"
erf.Ck              = 0.1
erf.Ce              = 0.7
erf.Ce_wall         = 3.9
erf.sigma_k         = 0.5

erf.theta_ref = 300.0

erf.rayleigh_damp_U = true
erf.rayleigh_damp_V = true
erf.rayleigh_damp_W = true
erf.rayleigh_damp_T = true
erf.rayleigh_zdamp  = 50.0
erf.rayleigh_dampcoef = 0.2

# INITIALIZATION
erf.init_type           = "input_sounding"
erf.sounding_type       = Ideal
erf.input_sounding_file = "stable_sounding.txt"

# INFLOW PERTURBATIONS
erf.fix_random_seed = 1
erf.perturbation_direction = 0 1 0 0 0 0
erf.perturbation_layers    = 3
erf.perturbation_offset    = 5
erf.perturbation_box_dims  = 3 3 3
erf.perturbation_klo       = 1
erf.perturbation_khi       = 29
erf.perturbation_Ug        = 13.5
erf.perturbation_nondimensional = 0.001
erf.perturbation_T_infinity     = 300
