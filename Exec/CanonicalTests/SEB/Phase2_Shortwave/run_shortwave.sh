#!/bin/bash
# Phase 2 regtest: shortwave with ray-cast shadowing.
#
#   ./run_shortwave.sh /path/to/erf_exec        # NP=4 by default
#
# Three decks on NP ranks: a fixed sun 30 deg and 70 deg above the western
# horizon, checked face by face against the analytic incidence and shadow
# height, and the solar mode at Boulder's solstice noon checked against the
# solar formulas. The zen60 deck is also run on one rank and its face dump
# must agree with the NP-rank one.
set -u
EXE=${1:?usage: run_shortwave.sh /path/to/erf_exec}
NP=${NP:-4}
rm -f ibseb_*.csv faces_*.csv; rm -rf plt0* chk0*
status=0
for v in zen60 zen20 solar; do
    echo "== $v ($NP ranks)"
    mpirun -np $NP "$EXE" inputs_$v > run_$v.log 2>&1 || { echo "run failed (see run_$v.log)"; exit 1; }
    grep "\[IBSEB DEBUG\] lev=0 sun:" run_$v.log | tail -1 | cut -c1-150
    case $v in
        zen60) python3 check_shortwave.py faces_zen60 60 270 800 100 0.3 0.2 || status=1 ;;
        zen20) python3 check_shortwave.py faces_zen20 20 270 800 100 0.3 0.2 || status=1 ;;
        solar) python3 check_shortwave.py faces_solar solar run_solar.log || status=1 ;;
    esac
done
echo "== zen60 on 1 rank against $NP ranks"
mv faces_zen60.rank0.csv faces_zen60_np.csv 2>/dev/null; cat faces_zen60.rank*.csv 2>/dev/null | grep -v "^i," >> faces_zen60_np.csv
"$EXE" inputs_zen60 > run_zen60_np1.log 2>&1 || { echo "run failed"; exit 1; }
python3 - <<'PY' || status=1
import numpy as np
def load(fn):
    a = np.loadtxt(fn, delimiter=",", skiprows=1); return a[np.lexsort((a[:,4], a[:,3], a[:,2], a[:,1], a[:,0]))]
a = load("faces_zen60_np.csv"); b = load("faces_zen60.rank0.csv")
ok = a.shape == b.shape and np.allclose(a, b, rtol=0, atol=1e-9)
print(f"  rank independence of the face dump: {'PASS' if ok else 'FAIL'} ({a.shape[0]} vs {b.shape[0]} faces)")
raise SystemExit(0 if ok else 1)
PY
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
