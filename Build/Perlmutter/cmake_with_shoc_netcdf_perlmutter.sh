#!/bin/bash

# Load the needed modules
module load gcc-native cmake cray-mpich cray-libsci cray-hdf5-parallel cray-netcdf-hdf5parallel

# Deactivate GPU aware MPI for CPU build
export MPICH_GPU_SUPPORT_ENABLED=0
export CRAY_ACCEL_TARGET=none

# Deduce the lib paths and files with $(CC/cc/ftn --cray-print-opts=libs)
CRAY_LIBS_CLEAN=$(CC --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')
CRAY_LIBS_CLEAN="$CRAY_LIBS_CLEAN $(cc --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')"
CRAY_LIBS_CLEAN="$CRAY_LIBS_CLEAN $(ftn --cray-print-opts=libs | sed 's/-Wl,--as-needed,//g; s/,--no-as-needed//g; s/,-l/ -l/g')"

# Configure and build
cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_C_COMPILER=cc \
      -DCMAKE_CXX_COMPILER=CC \
      -DCMAKE_C_FLAGS="$(cc --cray-print-opts=cflags)" \
      -DCMAKE_CXX_FLAGS="$(CC --cray-print-opts=cflags)" \
      -DCMAKE_CUDA_FLAGS="$(CC --cray-print-opts=cflags)" \
      -DCMAKE_CXX_STANDARD_LIBRARIES="-lmpi_gnu_123" \
      -DCMAKE_CUDA_STANDARD_LIBRARIES="-lmpi_gnu_123" \
      -DCMAKE_EXE_LINKER_FLAGS="-Wl,--no-as-needed $CRAY_LIBS_CLEAN" \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_TESTS:BOOL=ON \
      -DERF_ENABLE_SHOC:BOOL=ON \
      -DERF_ENABLE_HDF5:BOOL=ON \
      -DERF_ENABLE_NETCDF:BOOL=ON \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      .. && make -j8
