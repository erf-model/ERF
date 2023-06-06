#!/bin/bash

module load PrgEnv-gnu craype-accel-amd-gfx90a cray-mpich rocm cmake ccache ninja git
which clang
which clang++
which hipcc
cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_HIP:BOOL=ON \
      -DAMReX_AMD_ARCH=gfx90a \
      -DERF_ENABLE_TESTS:BOOL=ON \
      -DERF_ENABLE_ALL_WARNINGS:BOOL=ON \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DCMAKE_C_COMPILER=$(which hipcc) \
      -DCMAKE_CXX_COMPILER=$(which hipcc) \
      -DCMAKE_CXX_STANDARD=17 \
      .. && make -j8
