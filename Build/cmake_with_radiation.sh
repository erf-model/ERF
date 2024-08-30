#!/bin/bash

# Example CMake config script for an OSX laptop with OpenMPI

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_TESTS:BOOL=ON \
      -DERF_ENABLE_RRTMGP:BOOL=ON \
      -DERF_ENABLE_NETCDF:BOOL=ON \
      -DERF_ENABLE_HDF5:BOOL=ON \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      .. && make -j8
