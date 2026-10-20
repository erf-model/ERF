# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# A small dry-ish convective box over the simple land model, sized so that a
# restart test of SLM costs a second rather than the hours the SLM_AWAKEN and
# SLM_CASS_SAMRadiation decks cost -- both of those are gated behind
# ERF_TEST_ENABLE_EXTRA_LSM_TESTS and need input files this repository does not
# carry, so nothing exercised an SLM restart at all before this deck.
#
# What it is for: SLM keeps prognostic state that ERF does not expose through
# lsm_data -- the canopy, skin, ground-skin and canopy-air-space temperatures,
# the canopy water store -- and that state is checkpointed. erf.plot_lsm writes
# it to plt_lsm_2D_*, which is what the restart test compares; without that the
# comparison would see only the atmosphere and the soil columns and pass while
# the canopy state restarted from its initialization value. That is the failure
# originally reported on issue 4225 (SLM stopping on t_canop > tfriz after a
# restart that re-made the level-0 grids).
#
# Four soil layers rather than the nine the CASS case uses, and 8 x 8 x 16
# cells: enough for the surface energy balance to evolve every one of those
# fields over 20 steps, which is what makes the comparison mean anything, and
# small enough to run in CI.
#
erf.prob_name = "ABL"

max_step  = 20
stop_time = 1.0e6

amrex.fpe_trap_invalid = 1
fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =   0.   0.    0.
geometry.prob_hi     = 400. 400.  400.
amr.n_cell           =   8    8    16
geometry.is_periodic = 1 1 0

amr.max_level = 0

# The land model needs a surface layer to exchange with; Moeng fluxes rather
# than prescribed ones so the exchange actually responds to the canopy state.
zlo.type = "surface_layer"
zhi.type = "SlipWall"

erf.surface_layer.flux_type = "moeng"
erf.most.z0   = 0.1
erf.most.zref = 12.5

# SLM requires a moisture model
erf.moisture_model     = "Kessler"
erf.land_surface_model = "SLM"

slm.nsoil     = 4
slm.soil_dz   = 0.05 0.10 0.30 0.55
slm.landtype0 = 10      # evergreen forest
slm.LAI0      = 2.0
slm.clay0     = 13.0
slm.sand0     = 17.0
slm.sw0       = 0.60 0.65 0.70 0.80
slm.st0       = 300.0 299.8 299.5 299.0
slm.relax_hgt = 0.0 0.0 0.0 1.0
slm.soiltnudging = true
slm.soilwnudging = true
slm.tausoil      = 86400.0
slm.tabs_s   = 0.0
slm.t00      = 300.0
slm.z0_soil  = 0.0387
slm.mws_mx0  = 50.0
slm.Rc_max   = 5000.0
slm.T_opt    = 298.0

erf.use_gravity = true
erf.fixed_dt    = 0.5

erf.v = 1
amr.v = 1

erf.check_int   = -1
erf.plot_int_1  = -1
erf.plot_file_1 = plt
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta

# plt_lsm_* alongside plotfile 1: this is how the land state, including the
# fields lsm_data does not carry, reaches a comparison.
erf.plot_lsm = true

erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"

erf.les_type        = "Smagorinsky"
erf.Cs              = 0.1
erf.molec_diff_type = "None"
