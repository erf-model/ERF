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

cmake -DCMAKE_INSTALL_PREFIX:PATH=$ERF_INSTALL_DIR \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_SYCL:BOOL=ON \
      -DERF_ENABLE_TESTS:BOOL=ON \
      -DERF_ENABLE_ALL_WARNINGS:BOOL=ON \
      -DERF_ENABLE_RRTMGP:BOOL=ON \
      -DERF_ENABLE_NETCDF:BOOL=ON \
      -DERF_ENABLE_HDF5:BOOL=ON \
      -DERF_ENABLE_FFTW:BOOL=ON \
      -DERF_ENABLE_MASA:BOOL=OFF \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      --log-context --log-level STATUS \
      -B $ERF_BUILD_DIR -S $ERF_SOURCE_DIR && \
cmake --build $ERF_BUILD_DIR -j10 && \
cmake --install $ERF_BUILD_DIR --prefix=$ERF_INSTALL_DIR
