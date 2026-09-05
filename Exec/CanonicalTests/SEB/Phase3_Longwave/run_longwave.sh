#!/bin/bash
# Phase 3 regtest: view fractions and longwave.
#
#   ./run_longwave.sh /path/to/erf_exec        # NP=4 by default
#
# Two decks on NP ranks (sky longwave fixed, or gray from the air
# temperature) checked face by face against an independent hemisphere
# sampling and the longwave formulas; the fixed deck also on one rank, whose
# dump must equal the NP-rank one.
set -u
EXE=${1:?usage: run_longwave.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_*.csv faces_*.csv; rm -rf plt0* chk0*
status=0
for v in fixed gray; do
    echo "== $v ($NP ranks)"
    mpirun -np $NP "$EXE" inputs_$v > run_$v.log 2>&1 || { echo "run failed (see run_$v.log)"; exit 1; }
    grep "view fractions" run_$v.log | tail -1 | cut -c1-140
    case $v in
        fixed) python3 check_longwave.py faces_fixed fixed 300 0.9 0.95 300 16 8 || status=1 ;;
        gray)  python3 check_longwave.py faces_gray  gray  0.83 0.9 0.95 300 16 8 || status=1 ;;
    esac
done
echo "== fixed on 1 rank against $NP ranks"
cat faces_fixed.rank*.csv | awk 'NR==1 || !/^i,/' > faces_fixed_np.csv
"$EXE" inputs_fixed > run_fixed_np1.log 2>&1 || { echo "run failed"; exit 1; }
python3 - <<'PY' || status=1
import numpy as np
def load(fn):
    a = np.loadtxt(fn, delimiter=",", skiprows=1); return a[np.lexsort((a[:,4], a[:,3], a[:,2], a[:,1], a[:,0]))]
a = load("faces_fixed_np.csv"); b = load("faces_fixed.rank0.csv")
ok = a.shape == b.shape and np.allclose(a, b, rtol=0, atol=1e-9)
print(f"  rank independence of the face dump: {'PASS' if ok else 'FAIL'} ({a.shape[0]} vs {b.shape[0]} faces)")
raise SystemExit(0 if ok else 1)
PY
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
