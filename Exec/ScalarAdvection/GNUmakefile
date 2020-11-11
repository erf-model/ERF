AMREX_HOME ?= $(HOME)/Codes/ERF/Submodules/AMReX

# AMReX
DIM = 3
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_DPCPP = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# ERF
Eos_dir := GammaLaw
Reactions_dir := Null
Transport_dir := Constant

# GNU Make
Bpack := ./Make.package
Blocs := .
ERF_HOME := ../..
include $(ERF_HOME)/Exec/Make.ERF
