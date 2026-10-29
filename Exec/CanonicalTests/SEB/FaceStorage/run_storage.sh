#!/bin/bash
# Regression test: face storage.
#
#   ./run_storage.sh /path/to/erf_exec
#
# Runs the deck on one rank and on four ranks, checks the face counts against
# the blanking mask of the plotfile and that they agree between rank counts,
# then a checkpoint / restart round trip of the prognostic balance whose CSV
# rows and per-face dumps must match the straight run's, on the checkpoint's
# rank count and on one rank; last the two mismatch aborts.
set -u
EXE=${1:?usage: run_storage.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_*.csv faces_*.csv; rm -rf plt0* chk0*

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

echo "== straight prognostic run to step 4 ($NP ranks): the reference of the restart legs"
rm -rf plt0* chk0*
mpirun -np $NP "$EXE" inputs_straight_prognostic > run_straight_prognostic.log 2>&1 || { echo "run failed (see run_straight_prognostic.log)"; exit 1; }
grep "\[IBSEB\] lev=0" run_straight_prognostic.log | tail -1
echo "== checkpoint at step 2, restart to step 4 ($NP ranks)"
rm -rf plt0* chk0*
mpirun -np $NP "$EXE" inputs_chk > run_chk.log 2>&1 || { echo "chk run failed"; exit 1; }
[ -f chk00002/Level_0/IBSEBFaces ] || { echo "no IBSEBFaces in checkpoint: FAIL"; exit 1; }
# One record per face: the header line names the face count, then an 8-byte
# marker word, then 8 bytes of key and 8 x (2 + layers) bytes of values per
# face, nothing else.
python3 - <<'PY' || exit 1
import os
hdr = open("chk00002/Level_0/IBSEBFaces", "rb").readline()
tag, ver, nface, nl, nx, ny, nz = hdr.split()
expect = len(hdr) + 8 + int(nface) * 8 * (3 + int(nl))
size = os.path.getsize("chk00002/Level_0/IBSEBFaces")
ok = tag == b"IBSEBFaces" and nface == b"2056" and nl == b"4" and size == expect
print(f"IBSEBFaces: {nface.decode()} faces x (2 + {nl.decode()}) values, {size} bytes" + (" (the records alone): PASS" if ok else f" but expected {expect}: FAIL"))
raise SystemExit(0 if ok else 1)
PY
mpirun -np $NP "$EXE" inputs_restart > run_restart.log 2>&1 || { echo "restart run failed"; exit 1; }
grep -q "Face state restored" run_restart.log && echo "state restored: yes" || { echo "state restored: FAIL"; exit 1; }
# Every column of the row must match exactly, including the two read from the
# atmosphere (net longwave through the air temperature, the sensible flux
# through the wind), which a restart of the immersed forcing once moved.
python3 - <<'PY' || exit 1
s = open("ibseb_straight_prognostic.csv").read().strip().splitlines(); r = open("ibseb_restart.csv").read().strip().splitlines()
hdr = s[0].split(","); a = s[-1].split(","); b = r[-1].split(",")
ok = True
for h, x, y in zip(hdr, a, b):
    if x != y: ok = False; print(f" {h}: {x} vs {y} (must be exact)")
print("restart CSV row matches straight (exact): " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
PY
# And every face: skin temperature, slab layers and fluxes of the restart's step 4
# against the straight run's, record by record.
python3 compare_faces.py faces_straight_prognostic faces_restart || { echo "restart faces match straight: FAIL"; exit 1; }
echo "restart faces match straight (every column of every face): PASS"
echo "== the same restart on 1 rank (every face finds its record in the $NP-rank checkpoint)"
"$EXE" inputs_restart erf.ibseb.csv_file=ibseb_restart_np1.csv erf.ibseb.dump_faces_file=faces_restart_np1 > run_restart_np1.log 2>&1 || { echo "restart run failed"; exit 1; }
# Steps 3 and 4 run on another rank count, so the atmosphere differs in its last
# digit (ERF's own solve; the straight run on 1 rank differs from the one on $NP
# ranks the same way) and every column read from it follows. The comparison is
# therefore to 1e-8 relative, with the balance residual (|balance| after the
# Newton solve, ~1e-8 of terms of hundreds of W/m2) within 1e-6 W/m2; a lost or
# misplaced record would move the skin temperature by kelvins.
python3 - <<'PY' || exit 1
s = open("ibseb_restart.csv").read().strip().splitlines(); r = open("ibseb_restart_np1.csv").read().strip().splitlines()
hdr = s[0].split(","); a = s[-1].split(","); b = r[-1].split(",")
ok = True
for h, x, y in zip(hdr, a, b):
    if x == y: continue
    fx, fy = float(x), float(y)
    tol = 1e-6 if h.startswith("resid_max") else 1e-8 * max(abs(fx), abs(fy))
    if abs(fx - fy) > tol: ok = False; print(f" {h}: {x} vs {y}")
print("restart on 1 rank equals the restart on the checkpoint's rank count (to 1e-8, resid_max within 1e-6): " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
PY
python3 compare_faces.py faces_restart faces_restart_np1 --rtol 1e-8 --resid-atol 1e-6 || { echo "restart faces on 1 rank match: FAIL"; exit 1; }
echo "restart faces on 1 rank match the $NP-rank restart (every face, to 1e-8): PASS"
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
# The rotation keeps the face count (2056), so the count check passes and the
# per-face lookup is what must stop the run.
if grep -q "faces of this blanking have no record there" run_restart_rotated.log; then echo "checkpoint of other buildings rejected: PASS"; else echo "checkpoint of other buildings rejected: FAIL (see run_restart_rotated.log)"; exit 1; fi
echo "ALL PASS"
