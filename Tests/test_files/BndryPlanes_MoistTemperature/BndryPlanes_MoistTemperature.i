# Reads the planes BndryPlanes_MoistTemperature_Write.i writes (temperature, density, velocity
# and qv) at the x faces of its box, starting from the same uniform state.  Inflow and interior
# then agree exactly, so theta must stay 300 K to round-off.  The temperature planes hold the
# moist temperature, and the reader must convert them back to theta with the qv planes: without
# qv the inflow theta is 0.69 K too warm (erf-model/ERF#4290).
erf.prob_name = "ABL"

max_step = 50

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo =  256.  256.     0.
geometry.prob_hi =  768.  768.  1024.
amr.n_cell       =    8     8     16

geometry.is_periodic = 0 1 0

xlo.type = "Inflow"
xhi.type = "Outflow"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.substepping_type = None
erf.fixed_dt         = 2.0e-2

erf.v              = 1
amr.v              = 1

amr.max_level       = 0

erf.check_int       = -1
erf.plot_file_1     = plt
erf.plot_int_1      = 50
erf.plot_vars_1     = density x_velocity y_velocity z_velocity theta qv

erf.use_gravity = true

erf.molec_diff_type = "None"
erf.les_type        = "None"

erf.moisture_model = "SatAdj"

erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"

erf.input_bndry_planes    = 1
erf.bndry_file            = "BndryFiles"
erf.bndry_input_var_names = temperature density velocity qv
