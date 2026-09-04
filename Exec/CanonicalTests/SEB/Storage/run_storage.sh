#!/bin/bash
# Phase 1 regtest: face storage.
#
#   ./run_storage.sh /path/to/erf_exec
#
# Runs the deck on one rank and on four ranks, checks the face counts against
# the blanking mask of the plotfile and that they agree between rank counts,
# then a checkpoint / restart round trip whose CSV rows must match the
# straight run's.
set -u
EXE=${1:?usage: run_storage.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_*.csv; rm -rf plt0* chk0*

echo "== straight, 1 rank"
"$EXE" inputs_straight > run_straight_np1.log 2>&1 || { echo "run failed (see run_straight_np1.log)"; exit 1; }
grep "\[IBSEB\] lev=0" run_straight_np1.log | tail -1
python3 check_storage.py plt00004 run_straight_np1.log 300.0 || exit 1
line1=$(grep "\[IBSEB\] lev=0" run_straight_np1.log | tail -1 | sed 's/.*faces=/faces=/; s/ T_skin.*//')
cp ibseb_straight.csv ibseb_straight_np1.csv

echo "== straight, $NP ranks"
rm -f ibseb_straight.csv; rm -rf plt0*
mpirun -np $NP "$EXE" inputs_straight > run_straight_np${NP}.log 2>&1 || { echo "run failed (see run_straight_np${NP}.log)"; exit 1; }
lineN=$(grep "\[IBSEB\] lev=0" run_straight_np${NP}.log | tail -1 | sed 's/.*faces=/faces=/; s/ T_skin.*//')
python3 check_storage.py plt00004 run_straight_np${NP}.log 300.0 || exit 1
[ "$line1" = "$lineN" ] && echo "rank independence: PASS ($line1)" || { echo "rank independence: FAIL"; echo " 1: $line1"; echo " $NP: $lineN"; exit 1; }

echo "== checkpoint at step 2, restart to step 4 ($NP ranks)"
rm -rf plt0* chk0*
mpirun -np $NP "$EXE" inputs_chk > run_chk.log 2>&1 || { echo "chk run failed"; exit 1; }
grep -q "IBSEBState" chk00002/Level_0/IBSEBState_H 2>/dev/null || ls chk00002/Level_0 | grep -q IBSEBState || { echo "no IBSEBState in checkpoint: FAIL"; exit 1; }
mpirun -np $NP "$EXE" inputs_restart > run_restart.log 2>&1 || { echo "restart run failed"; exit 1; }
grep -q "Face state restored" run_restart.log && echo "state restored: yes" || { echo "state restored: FAIL"; exit 1; }
s=$(tail -1 ibseb_straight.csv); r=$(tail -1 ibseb_restart.csv)
[ "$s" = "$r" ] && echo "restart CSV row matches straight: PASS" || { echo "restart CSV row: FAIL"; echo " straight: $s"; echo " restart:  $r"; exit 1; }
echo "ALL PASS"
