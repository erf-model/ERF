#!/bin/bash
# Phase 5 regtest: slab conduction and materials.
#
#   ./run_ground.sh /path/to/erf_exec        # NP=4 by default
set -u
EXE=${1:?usage: run_ground.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_*.csv faces_*.csv; rm -rf plt0* chk0*
status=0
echo "== thick slab, semi-infinite response ($NP ranks, 50 s)"
mpirun -np $NP "$EXE" inputs_thick > run_thick.log 2>&1 || { echo "run failed (see run_thick.log)"; exit 1; }
python3 check_ground.py thick run_thick.log 1.5 2.0e6 20 0.25 faces_thick 0.2 200 50.0 || status=1
echo "== thin slab, steady linear profile ($NP ranks, 100 s)"
mpirun -np $NP "$EXE" inputs_thin > run_thin.log 2>&1 || { echo "run failed (see run_thin.log)"; exit 1; }
python3 check_ground.py thin faces_thin 1.0 0.02 20 || status=1
echo "== materials by building ($NP ranks)"
mpirun -np $NP "$EXE" inputs_materials > run_materials.log 2>&1 || { echo "run failed (see run_materials.log)"; exit 1; }
grep "IBSEB DEBUG\]   material" run_materials.log | cut -c1-120
python3 check_ground.py materials faces_materials materials.csv || status=1
echo "== slab through a checkpoint ($NP ranks)"
rm -rf chk0*
mpirun -np $NP "$EXE" inputs_chk > run_chk.log 2>&1 || { echo "chk run failed"; exit 1; }
mpirun -np $NP "$EXE" inputs_restart > run_restart.log 2>&1 || { echo "restart run failed"; exit 1; }
python3 check_ground.py restart faces_thin faces_restart || status=1
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
