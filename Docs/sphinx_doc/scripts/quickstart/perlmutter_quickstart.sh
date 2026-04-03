#!/usr/bin/env bash
set -euo pipefail

git clone --recursive https://github.com/erf-model/ERF.git
cd ERF
source Build/machines/perlmutter_erf.profile

mkdir -p build
cd build
../Build/cmake_with_kokkos_many_cuda.sh
make -j

cd Exec
srun -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most
