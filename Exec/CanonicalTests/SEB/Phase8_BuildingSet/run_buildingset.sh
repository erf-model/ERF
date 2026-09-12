#!/bin/bash
# Canonical building-set morning: run, check, plot.
#
#   ./run_buildingset.sh /path/to/erf_exec        # NP=4 by default, about an hour
set -u
EXE=${1:?usage: run_buildingset.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_set.csv; rm -rf faces plt* chk*; mkdir -p faces
echo "== building set, 6 h from 05:00 ($NP ranks)"
mpirun -np $NP "$EXE" inputs > run_set.log 2>&1 || { echo "run failed (see run_set.log)"; exit 1; }
python3 check_buildingset.py ibseb_set.csv faces/set run_set.log && echo "ALL PASS" || { echo "SOME CHECKS FAILED"; exit 1; }
