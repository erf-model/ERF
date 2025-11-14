#!/bin/bash

#Example cmake configuration script that assumes cray detection

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install_erf \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_FFT:BOOL=ON \
      -DERF_ENABLE_NETCDF:BOOL=ON \
      -DERF_ENABLE_HDF5:BOOL=ON \
      -DERF_ENABLE_RRTMGP:BOOL=ON \
      -DERF_ENABLE_SHOC:BOOL=OFF \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_CUDA:BOOL=ON \
      -DERF_ENABLE_HIP:BOOL=OFF \
      -DERF_ENABLE_SYCL:BOOL=OFF \
      -DERF_ENABLE_TESTS:BOOL=ON \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      --log-context \
      -B build_erf ..

cmake --build build_erf -j10 -v
cmake --install build_erf --prefix=install_erf
