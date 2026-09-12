#!/bin/bash
# Phase 6 regtest: the prognostic skin temperature.
#
#   ./run_prognostic.sh /path/to/erf_exec        # NP=4 by default
set -u
EXE=${1:?usage: run_prognostic.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_*.csv faces_*.csv; rm -rf plt0* chk0*
status=0
echo "== closure under a fixed sun ($NP ranks, 200 steps)"
mpirun -np $NP "$EXE" inputs_closure > run_closure.log 2>&1 || { echo "run failed (see run_closure.log)"; exit 1; }
grep "\[IBSEB\] lev=0 step=199" run_closure.log | cut -c1-400
python3 check_prognostic.py closure faces_closure 0.5 293.0 || status=1
echo "== external flux with the bound raised ($NP ranks, 100 steps)"
mpirun -np $NP "$EXE" inputs_qext > run_qext.log 2>&1 || { echo "run failed (see run_qext.log)"; exit 1; }
grep "\[IBSEB\] lev=0 step=99" run_qext.log | cut -c1-400
python3 check_prognostic.py qext faces_qext 0.5 293.0 3000.0 380.0 || status=1
echo "== through a checkpoint ($NP ranks)"
mpirun -np $NP "$EXE" inputs_chk > run_chk.log 2>&1 || { echo "chk run failed"; exit 1; }
mpirun -np $NP "$EXE" inputs_restart > run_restart.log 2>&1 || { echo "restart run failed"; exit 1; }
python3 check_prognostic.py restart faces_closure.step000199 faces_restart.step000199 || status=1
echo "== sunrise over the cube ($NP ranks, 10800 steps)"
mpirun -np $NP "$EXE" inputs_solar > run_solar.log 2>&1 || { echo "run failed (see run_solar.log)"; exit 1; }
python3 check_prognostic.py solar faces_solar || status=1
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
