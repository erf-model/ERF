# Reads the planes of BndryPlanes_Ingest_Writer.i on its x-low inflow face. The
# run starts at rest and dry, so the plotfile shows whether a field was taken
# from the planes; Tests/RunBndryPlanesIngest.cmake sets erf.bndry_input_var_names.
erf.prob_name = "ABL"
max_step = 20
amrex.fpe_trap_invalid = 1
geometry.prob_lo =  256.  256.     0.
geometry.prob_hi =  768.  768.  1024.
amr.n_cell       =    8     8     16
amr.max_grid_size   = 4    # four boxes, so the ingested faces cross box seams
amr.blocking_factor = 4
geometry.is_periodic = 0 1 0
xlo.type = "Inflow"
xhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"
erf.substepping_type = None
erf.fixed_dt         = 0.05
amr.max_level       = 0
erf.check_int       = -1
erf.plot_int_1      = 20
erf.plot_vars_1     = density x_velocity theta qv
erf.use_gravity = true
erf.molec_diff_type = "None"
erf.init_type = "input_sounding"
erf.input_sounding_file = "sounding_dry"
erf.moisture_model = "Kessler_NoRain"
erf.input_bndry_planes = 1
erf.bndry_file = "BndryFiles"
erf.bndry_input_var_names = velocity  # set by the test
