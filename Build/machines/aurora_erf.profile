#!/bin/bash
# Aurora (ALCF) environment profile for ERF
# See also: Docs/sphinx_doc/aurora_build_run.rst

# Load required modules
module load mpich/opt/4.2.3-intel
module load hdf5/1.14.6
module load python/3.10.14
module load cmake

# Intel compilers for MPICH wrappers
export MPICH_CC=icx
export MPICH_CXX=icpx
export MPICH_FC=ifx
export MPICH_F90=ifx

# NetCDF: User must provide NETCDF_DIR for NetCDF-enabled builds.
# C-only NetCDF is sufficient for ERF NetCDF I/O.
# NetCDF-Fortran support is required only if enabling features that use it (e.g., Noah-MP).
# For an example parallel NetCDF build with Fortran support, see:
#   https://github.com/Lab-Notebooks/ERF-Noah-Coupling/tree/main/software/netcdf
#
# Quick checks:
#   ls $NETCDF_DIR/include/netcdf.h       # C headers
#   ls $NETCDF_DIR/include/netcdf.mod     # Fortran module (if Fortran support built)
#   ls $NETCDF_DIR/lib/libnetcdf*         # C library
#   ls $NETCDF_DIR/lib/libnetcdff*        # Fortran library (if Fortran support built)

# GPU-aware MPI (uncomment to enable)
# export MPICH_GPU_SUPPORT_ENABLED=1