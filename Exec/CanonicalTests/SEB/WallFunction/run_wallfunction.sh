#!/bin/bash
# Regression test: the wall function beyond neutral.
#
#   ./run_wallfunction.sh /path/to/erf_exec        # NP=4 by default, a few minutes
set -u
EXE=${1:?usage: run_wallfunction.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_*.csv faces_*.csv; rm -rf plt0* chk0*
status=0
for v in neutral deardorff stability bulkri; do
    echo "== $v ($NP ranks, 600 steps)"
    mpirun -np $NP "$EXE" inputs_$v > run_$v.log 2>&1 || { echo "run failed (see run_$v.log)"; exit 1; }
    grep "\[IBSEB\] lev=0 step=600" run_$v.log | sed 's/.*T_skin_min/T_skin_min/' | cut -c1-260
done
python3 check_wallfunction.py neutral   faces_neutral faces_deardorff || status=1
python3 check_wallfunction.py deardorff faces_deardorff 0.5 1000.0 1.2 || status=1
python3 check_wallfunction.py stability faces_stability faces_deardorff 1000.0 1.2 || status=1
python3 check_wallfunction.py bulkri    run_bulkri.log faces_bulkri 95.0 || status=1
echo "== deardorff through a checkpoint at step 300 ($NP ranks)"
mpirun -np $NP "$EXE" inputs_deardorff_chk > run_deardorff_chk.log 2>&1 || { echo "chk run failed"; exit 1; }
mpirun -np $NP "$EXE" inputs_deardorff_restart > run_deardorff_restart.log 2>&1 || { echo "restart run failed"; exit 1; }
python3 ../PrognosticSkin/check_prognostic.py restart faces_deardorff.step000600 faces_deardorff_restart.step000600 || status=1
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
