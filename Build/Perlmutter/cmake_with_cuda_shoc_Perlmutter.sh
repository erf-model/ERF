#!/bin/bash

module load gcc-native cmake cray-mpich cray-libsci cray-hdf5-parallel cray-netcdf-hdf5parallel

# NOTE: $(CC --cray-print-opts=libs) can be used to deduce libmpi_gnu_123.so

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_PREFIX_PATH:PATH=${CUDATOOLKIT_HOME}/../../ \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DCMAKE_C_COMPILER=cc \
      -DCMAKE_CXX_COMPILER=CC \
      -DCMAKE_C_FLAGS="$(cc --cray-print-opts=cflags)" \
      -DCMAKE_CXX_FLAGS="$(CC --cray-print-opts=cflags)" \
      -DCMAKE_CUDA_FLAGS="$(CC --cray-print-opts=cflags)" \
      -DCMAKE_CXX_STANDARD_LIBRARIES="-lmpi_gnu_123" \
      -DCMAKE_CUDA_STANDARD_LIBRARIES="-lmpi_gnu_123" \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_TESTS:BOOL=ON \
      -DERF_ENABLE_CUDA:BOOL=ON \
      -DERF_ENABLE_NVHPC:BOOL=ON \
      -DERF_ENABLE_EKAT:BOOL=ON \
      -DERF_ENABLE_SHOC:BOOL=ON \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      .. && make -j8
