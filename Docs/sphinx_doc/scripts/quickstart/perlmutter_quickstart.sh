#!/usr/bin/env bash
set -euo pipefail

git clone --recursive https://github.com/erf-model/ERF.git
cd ERF
source Build/machines/perlmutter_erf.profile

mkdir -p build
cd build
../Build/cmake_with_kokkos_many_cuda.sh
make -j

# Optional developer smoke test (requires -DERF_ENABLE_TESTS:BOOL=ON at configure time):
# ctest -L regression -VV -j 4

cd Exec
srun -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most
