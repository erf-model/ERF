# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# MassSource.i with the rate read from a table, mass_source_ramp.txt: a ramp from 0 to 1.0e4 kg/s
# over 5 s, held to 8 s and zero after it, 55000 kg in all (rate and rate_file cannot both be
# given, so this is a deck of its own).  The rest is as in MassSource.i:
#
# A localized mass source (erf.mass_source) in a closed box: slip walls on all six faces, so
# nothing enters or leaves.  The source injects a constant 1.0e4 kg/s
# from the cells within 100 m of the center and less than 100 m above the ground, with the
# tracer on.  The domain mass ERF prints ("MASS") must therefore grow by exactly 1.0e4 kg/s
# times the elapsed time, and the tracer mass ("RHO SCALAR"), zero at the start, must equal it
# (Tests/check_mass_source.py).  Tests/CTestList.cmake runs it flat, over a hill, on a
# stretched mesh, with the grids split in z, with a rate table (mass_source_ramp.txt), a fixed
# theta and no tracer, with refined patches that hold all, none, the column above or half of the
# source region, and anelastic, where ERF must stop.
#
max_step = 20

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =    0.    0.    0.
geometry.prob_hi     = 1000. 1000.  500.
amr.n_cell           =   32    32    16
amr.max_grid_size    =   16    16    16

geometry.is_periodic = 0 0 0

xlo.type = "SlipWall"
xhi.type = "SlipWall"
ylo.type = "SlipWall"
yhi.type = "SlipWall"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

erf.compressibility_type = FullyCompressible

# MASS SOURCE
erf.mass_source.center = 500.0 500.0
erf.mass_source.radius = 100.0
erf.mass_source.height = 100.0
erf.mass_source.rate_file = "mass_source_ramp.txt"
erf.mass_source.tracer = true

# INITIALIZATION: air at rest, 300 K
erf.init_type           = "input_sounding"
erf.input_sounding_file = "input_sounding"

# TIME STEP CONTROL: dt = 0.5 s, acoustic substeps of 0.05 s (c dt_fast / dx = 0.55)
erf.fixed_dt           = 0.5
erf.fixed_mri_dt_ratio = 10

# DIAGNOSTICS & VERBOSITY: the integrals at every step, with every digit
erf.sum_interval = 1
erf.v            = 1
amr.v            = 1

# REFINEMENT (off here; the AMR test refines a patch around the source)
amr.max_level = 0

# CHECKPOINT AND PLOT FILES
erf.check_int   = -1
erf.plot_int_1  = -1
erf.plot_vars_1 = density rhoadv_0

# SOLVER CHOICE
erf.use_gravity = true
erf.buoyancy_type = 1

erf.les_type        = "None"
erf.molec_diff_type = "None"
