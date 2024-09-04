#!/bin/bash

# Example CMake config script for an OSX laptop with OpenMPI

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_CUDA:BOOL=ON \
      -DERF_ENABLE_TESTS:BOOL=ON \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCUDAToolkit_ROOT=/usr/local/cuda-12.4/bin \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      .. && make -j8
