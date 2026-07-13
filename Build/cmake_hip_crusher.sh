#!/bin/bash
# Defaults (customize here or set ERF_BUILD_DIR/ERF_SOURCE_DIR/ERF_INSTALL_DIR in environment)
# If ERF_HOME is set, use it as base for absolute paths
if [ -n "$ERF_HOME" ]; then
  : ${ERF_BUILD_DIR:="$ERF_HOME/build"}
  : ${ERF_SOURCE_DIR:="$ERF_HOME"}
  : ${ERF_INSTALL_DIR:="$ERF_HOME/install"}
else
  : ${ERF_BUILD_DIR:="."}
  : ${ERF_SOURCE_DIR:=".."}
  : ${ERF_INSTALL_DIR:="install"}
fi

echo "Source: $ERF_SOURCE_DIR | Build: $ERF_BUILD_DIR | Install: $ERF_INSTALL_DIR | PWD: $(pwd)"
echo "Customize: export ERF_BUILD_DIR=... ERF_SOURCE_DIR=... ERF_INSTALL_DIR=... or ERF_HOME=..."

module load PrgEnv-gnu craype-accel-amd-gfx90a cray-mpich rocm cmake ccache ninja git
which clang
which clang++
which hipcc
cmake -DCMAKE_INSTALL_PREFIX:PATH=$ERF_INSTALL_DIR \
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
      --log-context --log-level STATUS \
      -B $ERF_BUILD_DIR -S $ERF_SOURCE_DIR && \
cmake --build $ERF_BUILD_DIR -j10 && \
cmake --install $ERF_BUILD_DIR --prefix=$ERF_INSTALL_DIR
