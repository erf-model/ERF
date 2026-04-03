#!/usr/bin/env bash
set -euo pipefail

# Set NETCDF_DIR before running this script.
: "${NETCDF_DIR:?Set NETCDF_DIR to your NetCDF installation prefix}"

git clone --recursive https://github.com/erf-model/ERF.git
cd ERF
source Build/machines/aurora_erf.profile

mkdir -p build
cd build
../Build/cmake_with_kokkos_many_sycl.sh
make -j

cd Exec
mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most
