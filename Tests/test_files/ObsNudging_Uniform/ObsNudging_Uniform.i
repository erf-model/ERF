# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Observation nudging with a closed-form answer.
#
# The flow is uniform, with no gravity, no diffusion and no surface drag, and
# the one station has a horizontal radius far larger than the domain and a
# profile that spans its whole height, so every cell has weight 1 (to within
# 1e-7) and relaxes as the same ODE,
#
#     d(phi)/dt = -(phi - (a + b t)) / tau,
#
# towards a target that the file makes linear in time.  Its solution is
#
#     phi(t) = a + b t - b tau + (phi_0 - a + b tau) exp(-t / tau),
#
# and StationSeriesCheck compares the station series written below with it.
# The same deck with a sigma band checks that the target is the near edge of
# the band, mean - alpha sigma, while the value is below it.
#
erf.prob_name = "ABL"

max_step = 40

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =   0.   0.   0.
geometry.prob_hi     = 800. 800. 400.
amr.n_cell           =  16   16   16

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
# The acoustic substep resolves c dtau / dx = 347 * (1/12) / 50 = 0.58
erf.fixed_dt       = 1.0
erf.fixed_mri_dt_ratio = 12

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = -1
erf.v              = 1
amr.v              = 1

# REFINEMENT / REGRIDDING
amr.max_level      = 0

# CHECKPOINT FILES
erf.check_file     = chk
erf.check_int      = -1

# PLOTFILES
erf.plot_file_1    = plt
erf.plot_int_1     = 40
erf.plot_vars_1    = density x_velocity y_velocity z_velocity theta

# SOLVER CHOICE
erf.use_gravity     = false
erf.use_coriolis    = false
erf.molec_diff_type = "None"
erf.les_type        = "None"
erf.pbl_type        = "None"
erf.abl_driver_type = "None"

# INITIAL CONDITIONS
erf.init_type = Uniform
prob.rho_0 = 1.0
prob.A_0   = 1.0
prob.T_0   = 300.0
prob.U_0   = 4.0
prob.V_0   = 1.0
prob.W_0   = 0.0

# OBSERVATION NUDGING
erf.nudging_from_observations = true
erf.obs_nudging.stations          = mast
erf.obs_nudging.tau               = 20.0
erf.obs_nudging.horizontal_radius = 1.0e6
erf.obs_nudging.vertical_radius   = 25.0
erf.obs_nudging.sigma_factor      = 0.0
erf.obs_nudging.mast.file         = uniform_station.txt
erf.obs_nudging.mast.x            = 400.0
erf.obs_nudging.mast.y            = 400.0

# STATION OUTPUT, far from the nudging station
erf.station_names             = probe
erf.probe.field               = x_velocity y_velocity z_velocity theta
erf.probe.x                   = 125.0
erf.probe.y                   = 675.0
erf.probe.height_agl          = 37.5 212.5 362.5
erf.station_sampling_interval = 1
