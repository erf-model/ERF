#!/bin/bash
# Aurora (ALCF) environment profile for ERF
# See also: Docs/sphinx_doc/aurora_build_run.rst

# Load required modules
module load mpich/opt/4.2.3-intel
module load hdf5/1.14.6
module load netcdf-cxx4
module load netcdf-c
module load netcdf-fortran
module load python/3.10.14
module load cmake

# Intel compilers for MPICH wrappers
export MPICH_CC=icx
export MPICH_CXX=icpx
export MPICH_FC=ifx
export MPICH_F90=ifx

# NetCDF and HDF5 roots are set automatically by the modules above:
#   NETCDF_C_ROOT      — set by netcdf-c module
#   NETCDF_FORTRAN_ROOT — set by netcdf-fortran module
#   HDF5_ROOT           — set by hdf5 module
# If NETCDF_FORTRAN_ROOT is not set by the module, derive it:
if [[ -z "${NETCDF_FORTRAN_ROOT:-}" ]]; then
  export NETCDF_FORTRAN_ROOT=$(module show netcdf-fortran 2>&1 | \
    sed -n 's/.*setenv("NETCDF_FORTRAN_ROOT","\([^"]*\)").*/\1/p')
fi
#
# Quick checks:
#   ls $NETCDF_C_ROOT/include/netcdf.h        # C headers
#   ls $NETCDF_C_ROOT/lib64/libnetcdf*        # C library
#   ls $NETCDF_FORTRAN_ROOT/include/netcdf.mod # Fortran module
#   ls $NETCDF_FORTRAN_ROOT/lib64/libnetcdff* # Fortran library
#   ls $HDF5_ROOT/lib/libhdf5*                # HDF5 library

# GPU-aware MPI (uncomment to enable)
# export MPICH_GPU_SUPPORT_ENABLED=1
