#!/usr/bin/env bash
set -euo pipefail

git clone --recursive https://github.com/erf-model/ERF.git
cd ERF
source Build/machines/polaris_erf.profile

mkdir -p build
cd build
../Build/cmake_with_kokkos_many_cuda.sh
make -j

# Optional interactive allocation example:
# qsub -I -A <PROJECT> -q debug -l select=1 -l walltime=1:00:00
# Optional batch example:
# qsub submit_erf_polaris.pbs

cd Exec
mpiexec -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most
