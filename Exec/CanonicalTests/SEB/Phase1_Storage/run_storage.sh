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
# The row's geometry, skin and slab columns must match exactly; the columns
# read from the atmosphere (net longwave through the air temperature, the
# sensible flux through the wind) may differ by the immersed forcing's
# restart non-exactness of about 1e-5 relative (see Phase5_Ground).
python3 - <<'PY' || exit 1
s = open("ibseb_straight.csv").read().strip().splitlines(); r = open("ibseb_restart.csv").read().strip().splitlines()
hdr = s[0].split(","); a = s[-1].split(","); b = r[-1].split(",")
loose = {"LW_net_mean_Wm2", "H_mean_Wm2"}; ok = True
for h, x, y in zip(hdr, a, b):
    if h in loose:
        if abs(float(x) - float(y)) > 1e-3 * max(1.0, abs(float(x))): ok = False; print(f" {h}: {x} vs {y} (beyond 1e-3 relative)")
    elif x != y: ok = False; print(f" {h}: {x} vs {y} (must be exact)")
print("restart CSV row matches straight (exact, atmosphere columns to 1e-3): " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
PY
echo "ALL PASS"
