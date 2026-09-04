#!/bin/bash
# Immersed-boundary awareness of MRF and YSUNew.
#
#   ./run_pbl_ib.sh /path/to/erf_exec        # NP=4 by default, a few minutes
set -u
set -o pipefail
EXE=${1:?usage: run_pbl_ib.sh /path/to/erf_exec}
NP=${NP:-4}
rm -rf plt_* Backtrace.*
status=0
for d in mrf_off ysunew_off; do   # allowed to fail: that is the defect
    echo "== $d ($NP ranks, 200 steps; may fail, see the README)"
    mpirun -np $NP "$EXE" inputs_$d > run_$d.log 2>&1 || echo "   run with the switch off failed, as documented"
done
for d in mrf_on ysunew_on flat_off flat_on; do
    echo "== $d ($NP ranks, 200 steps)"
    mpirun -np $NP "$EXE" inputs_$d > run_$d.log 2>&1 || { echo "run failed (see run_$d.log)"; exit 1; }
done
python3 check_pbl_ib.py scheme plt_mrf_on00200 plt_mrf_off00200 run_mrf_off.log 2>&1 | grep -v "^yt" || status=1
python3 check_pbl_ib.py scheme plt_ysunew_on00200 plt_ysunew_off00200 run_ysunew_off.log 2>&1 | grep -v "^yt" || status=1
python3 check_pbl_ib.py flat plt_flat_off00200 plt_flat_on00200 2>&1 | grep -v "^yt" || status=1
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
