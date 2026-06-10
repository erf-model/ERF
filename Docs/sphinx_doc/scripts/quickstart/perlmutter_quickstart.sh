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
# Tip: use interactive allocations for debugging; use batch jobs for testing.
# Optional interactive allocation example:
# salloc -A <account_or_project_g> -N 1 -t 00:30:00
# Optional batch example:
# sbatch ../Docs/sphinx_doc/scripts/quickstart/run_perlmutter_erf.sbatch

cd Exec
srun -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most
