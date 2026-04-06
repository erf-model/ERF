#!/usr/bin/env bash
set -euo pipefail

git clone --recursive https://github.com/erf-model/ERF.git
cd ERF
source Build/machines/kestrel_erf.profile

mkdir -p build
cd build
cp ../Docs/sphinx_doc/scripts/quickstart/cmake_kestrel_ERF.sh .
chmod +x cmake_kestrel_ERF.sh
./cmake_kestrel_ERF.sh

# Optional developer smoke test (requires -DERF_ENABLE_TESTS:BOOL=ON at configure time):
# ctest -L regression -VV -j 4
# Tip: use interactive allocations for debugging; use batch jobs for testing.

# Optional interactive allocation example:
# salloc -p gpu-h100s --nodes=1 --gpus-per-node=4 --time=00:30:00
# Optional batch example:
# sbatch ../Docs/sphinx_doc/scripts/quickstart/run.erf.aw.job_arena

cd Exec
srun -n 4 ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most
