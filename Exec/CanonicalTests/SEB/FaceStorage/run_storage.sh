#!/bin/bash
# Regression test: face storage.
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
# restart non-exactness of about 1e-5 relative (see SlabConduction).
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
echo "== the same restart on 1 rank (the state field is redistributed from the $NP-rank checkpoint)"
"$EXE" inputs_restart erf.ibseb.csv_file=ibseb_restart_np1.csv > run_restart_np1.log 2>&1 || { echo "restart run failed"; exit 1; }
python3 - <<'PY' || exit 1
s = open("ibseb_restart.csv").read().strip().splitlines()[-1]; r = open("ibseb_restart_np1.csv").read().strip().splitlines()[-1]
ok = s == r
print("restart on 1 rank equals the restart on the checkpoint's rank count: " + ("PASS" if ok else "FAIL"))
if not ok: print(" " + s); print(" " + r)
raise SystemExit(0 if ok else 1)
PY
echo "== restart with erf.ibseb.n_slab_layers = 6 against the 4-layer checkpoint (must abort)"
mpirun -np $NP "$EXE" inputs_restart_layers > run_restart_layers.log 2>&1
if grep -q "the deck sets erf.ibseb.n_slab_layers = 6" run_restart_layers.log; then echo "mismatched slab layers rejected: PASS"; else echo "mismatched slab layers rejected: FAIL (see run_restart_layers.log)"; exit 1; fi
echo "== restart against a checkpoint of other buildings (must abort): the height map rotated by 16 rows"
python3 - <<'PY'
lines = open("skyscraper_5m_128x128.txt").read().rstrip("\n").split("\n")
n = 128 * 128; head, heights = lines[:-n], lines[-n:]
shift = 16 * 128
open("skyscraper_rotated.txt", "w").write("\n".join(head + heights[-shift:] + heights[:-shift]) + "\n")
PY
mpirun -np $NP "$EXE" inputs_restart erf.buildings_file_name=skyscraper_rotated.txt erf.ibseb.csv_file=ibseb_restart_rotated.csv > run_restart_rotated.log 2>&1
if grep -q "was written for a different building layout" run_restart_rotated.log; then echo "checkpoint of other buildings rejected: PASS"; else echo "checkpoint of other buildings rejected: FAIL (see run_restart_rotated.log)"; exit 1; fi
echo "ALL PASS"
