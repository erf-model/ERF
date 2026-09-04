#!/bin/bash
# Phase 8 regtest: the wall function beyond neutral.
#
#   ./run_wallfunction.sh /path/to/erf_exec        # NP=4 by default, a few minutes
set -u
EXE=${1:?usage: run_wallfunction.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_*.csv faces_*.csv; rm -rf plt0*
status=0
for v in neutral deardorff stability bulkri; do
    echo "== $v ($NP ranks, 600 steps)"
    mpirun -np $NP "$EXE" inputs_$v > run_$v.log 2>&1 || { echo "run failed (see run_$v.log)"; exit 1; }
    grep "\[IBSEB\] lev=0 step=599" run_$v.log | sed 's/.*T_skin_min/T_skin_min/' | cut -c1-260
done
python3 check_wallfunction.py neutral   faces_neutral faces_deardorff || status=1
python3 check_wallfunction.py deardorff faces_deardorff 0.5 1000.0 1.2 || status=1
python3 check_wallfunction.py stability faces_stability faces_deardorff 1000.0 1.2 || status=1
python3 check_wallfunction.py bulkri    run_bulkri.log faces_bulkri 95.0 || status=1
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
