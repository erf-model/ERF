#!/bin/bash
# Canonical isolated-building day: run, check, plot.
#
#   ./run_isolated.sh /path/to/erf_exec        # NP=4 by default, about two hours
set -u
EXE=${1:?usage: run_isolated.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_day.csv; rm -rf faces plt* chk*; mkdir -p faces
echo "== isolated building, 24 h ($NP ranks)"
mpirun -np $NP "$EXE" inputs > run_day.log 2>&1 || { echo "run failed (see run_day.log)"; exit 1; }
python3 check_isolated.py ibseb_day.csv faces/day 0.5 295.0 && echo "ALL PASS" || { echo "SOME CHECKS FAILED"; exit 1; }
