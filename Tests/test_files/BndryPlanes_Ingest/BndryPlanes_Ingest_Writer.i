# Writes boundary planes of temperature, velocity, density and qv on the box
# [256,768]^2 from a moist, periodic run (u = 10 m/s, qv = 5 g/kg) for
# BndryPlanes_Ingest_Reader.i (Tests/RunBndryPlanesIngest.cmake).
erf.prob_name = "ABL"
max_step = 24
amrex.fpe_trap_invalid = 1
geometry.prob_lo =    0.    0.     0.
geometry.prob_hi = 1024. 1024.  1024.
amr.n_cell       =   16    16     16
geometry.is_periodic = 1 1 0
zlo.type = "SlipWall"
zhi.type = "SlipWall"
erf.substepping_type = None
erf.fixed_dt         = 0.05
amr.max_level       = 0
erf.check_int       = -1
erf.plot_int_1      = -1
erf.use_gravity = true
erf.molec_diff_type = "None"
erf.init_type = "input_sounding"
erf.input_sounding_file = "sounding_moist"
erf.moisture_model = "Kessler_NoRain"
erf.output_bndry_planes = 1
erf.bndry_output_planes_interval = 2
erf.bndry_output_start_time = 0.0
erf.bndry_output_planes_file = "BndryFiles"
erf.bndry_output_var_names = temperature velocity density qv
erf.bndry_output_box_lo = 256. 256.
erf.bndry_output_box_hi = 768. 768.
