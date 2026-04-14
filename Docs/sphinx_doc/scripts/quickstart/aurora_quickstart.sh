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

# Optional developer smoke test (requires -DERF_ENABLE_TESTS:BOOL=ON at configure time):
# ctest -L regression -VV -j 4
# Tip: use interactive allocations for debugging; use batch jobs for testing.
# Optional interactive allocation example:
# qsub -I -A <PROJECT> -q debug -l select=1 -l walltime=1:00:00 -l filesystems=flare
# Optional batch example:
# qsub ../Docs/sphinx_doc/scripts/quickstart/submit_erf_aurora.pbs

cd Exec
mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most
