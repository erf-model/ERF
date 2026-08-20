# AMReX
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE       = FALSE
TINY_PROFILE  = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE   = FALSE
USE_GPROF     = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE

USE_CUDA = FALSE
USE_HIP  = FALSE
USE_SYCL = FALSE

# Enable debug flags (rich backtrace)
DEBUG = FALSE

TEST = TRUE
USE_ASSERTION = TRUE

# Enable NETCDF support (required for WRFInput or MetGrid initialization)
USE_NETCDF = TRUE

# Enable FFT Poisson solver
USE_FFT = FALSE

# Enable NOAH LSM
USE_NOAHMP = FALSE

# Enable particles (required for ParticleTests or MultiSpeciesBubble) 
USE_PARTICLES = FALSE

# GNU Make
Bpack := ./Make.package
Blocs := .

ERF_HOME := ../
ERF_PROBLEM_DIR = $(ERF_HOME)/Exec
include $(ERF_HOME)/Exec/Make.ERF
