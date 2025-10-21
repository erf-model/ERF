#!/bin/bash

# Set ERF_DIR relative to the script location
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ERF_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
echo "ERF_DIR set to: $ERF_DIR"

# Set EKAT_HOME
EKAT_HOME="${ERF_DIR}/Submodules/ekat"
echo "EKAT_HOME set to: $EKAT_HOME"

SRC_DIR=${ERF_DIR}/Submodules/ekat
cd ${SRC_DIR}

# Cannot build/install in the src dir
WORK_DIR=${SRC_DIR}/build
INSTALL_DIR=${WORK_DIR}/install
mkdir -p ${INSTALL_DIR}
cd ${WORK_DIR}

# Set Kokkos source dir
KOKKOS_SRC=${SRC_DIR}/extern/kokkos

# Run CMake

cmake \
  -DCMAKE_INSTALL_PREFIX="${EKAT_HOME}" \
  -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} \
  -DCMAKE_INSTALL_INCLUDEDIR="include/kokkos" \
  -DCMAKE_INSTALL_LIBDIR="lib" \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ENABLE_CUDA_LAMBDA=ON \
  -DKokkos_ARCH_AMPERE80=ON \
  -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_CXX_STANDARD_REQUIRED=ON \
  -DKokkos_ENABLE_CXX17=ON \
  -DCMAKE_CUDA_STANDARD=17 \
  -DCMAKE_CUDA_STANDARD_REQUIRED=ON \
  -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON \
  -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_CUDA_COMPILER=nvcc \
  "$KOKKOS_SRC"

# Build and install
make -j8
make install

# Return to ERF directory
cd "$ERF_DIR" || { echo "Failed to return to ERF directory"; exit 1; }
echo "Returned to ERF directory: $ERF_DIR"
echo "Kokkos build and installation complete at ${EKAT_HOME}/kokkos"



