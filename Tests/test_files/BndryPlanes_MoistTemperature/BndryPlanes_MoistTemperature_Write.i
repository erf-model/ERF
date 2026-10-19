# Writer run for BndryPlanes_MoistTemperature (Tests/RunFieldBounds.cmake runs it first, in the
# same folder).  A moist atmosphere at rest in the frame of a uniform 10 m/s wind -- theta 300 K
# and qv 5 g/kg at every height, in hydrostatic balance -- on a periodic domain, writing
# temperature, density, velocity and qv planes every step from the box [256,768]^2 into
# BndryFiles.  The temperature written is the moist temperature, computed with qv.
erf.prob_name = "ABL"

max_step = 60

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

geometry.prob_lo =    0.    0.     0.
geometry.prob_hi = 1024. 1024.  1024.
amr.n_cell       =   16    16     16

geometry.is_periodic = 1 1 0

zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.substepping_type = None
erf.fixed_dt         = 2.0e-2

erf.v              = 1
amr.v              = 1

amr.max_level       = 0

erf.check_int       = -1
erf.plot_int_1      = -1

erf.use_gravity = true

erf.molec_diff_type = "None"
erf.les_type        = "None"

erf.moisture_model = "SatAdj"

erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"

erf.output_bndry_planes = 1
erf.bndry_output_planes_interval = 1
erf.bndry_output_start_time = 0.0
erf.bndry_output_planes_file = "BndryFiles"
erf.bndry_output_var_names = temperature density velocity qv

erf.bndry_output_box_lo = 256. 256.
erf.bndry_output_box_hi = 768. 768.
