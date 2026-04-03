#!/usr/bin/env bash
set -euo pipefail

git clone --recursive https://github.com/erf-model/ERF.git
cd ERF
source Build/machines/frontier_erf.profile

mkdir -p build
cd build
../Build/cmake_with_kokkos_many_hip.sh
make -j

# Optional interactive allocation example:
# salloc -A <project> -p batch -N 1 -t 00:30:00
# Optional batch example:
# sbatch run_frontier_erf.sbatch

cd Exec
srun -n 8 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most
