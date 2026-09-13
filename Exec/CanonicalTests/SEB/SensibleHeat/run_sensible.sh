#!/bin/bash
# Phase 4 regtest: sensible heat through the wall function.
#
#   ./run_sensible.sh /path/to/erf_exec        # NP=4 by default
#
# inputs_diag: the flux is diagnosed but not applied; every face's H must
#   equal the wall-function formula recomputed from the dumped wind, density
#   and temperatures, on NP ranks and on one rank identically.
# inputs_couple: the flux is applied; the heat the faces added over 40 steps
#   must match the extra enthalpy of the air against the diagnostic run at
#   the same step (the two runs differ only in the applied flux).
# inputs_inflow: mass inflow / pressure outflow; must run, with the wake
#   warmer than the inflow (plot_sensible.py shows it).
set -u
EXE=${1:?usage: run_sensible.sh /path/to/erf_exec}
NP=${NP:-4}
CP=1004.5
rm -f ibseb_*.csv faces_*.csv; rm -rf plt0* plt_diag_* chk0*
status=0
echo "== diag ($NP ranks)"
mpirun -np $NP "$EXE" inputs_diag > run_diag.log 2>&1 || { echo "run failed (see run_diag.log)"; exit 1; }
grep "\[IBSEB\] lev=0" run_diag.log | tail -1 | sed 's/.*T_skin_min/T_skin_min/' | cut -c1-140
python3 check_sensible.py formula faces_diag 0.01 0.001 $CP || status=1
cat faces_diag.rank*.csv | awk 'NR==1 || !/^i,/' > faces_diag_np.csv
mv plt00040 plt_diag_00040
echo "== diag (1 rank)"
"$EXE" inputs_diag > run_diag_np1.log 2>&1 || { echo "run failed"; exit 1; }
python3 - <<'PY' || status=1
import numpy as np
def load(fn):
    a = np.loadtxt(fn, delimiter=",", skiprows=1); return a[np.lexsort((a[:,4], a[:,3], a[:,2], a[:,1], a[:,0]))]
a = load("faces_diag_np.csv"); b = load("faces_diag.rank0.csv")
ok = a.shape == b.shape and np.allclose(a, b, rtol=1e-12, atol=1e-9)
print(f"  rank independence of the face dump: {'PASS' if ok else 'FAIL'} ({a.shape[0]} vs {b.shape[0]} faces)")
raise SystemExit(0 if ok else 1)
PY
echo "== couple ($NP ranks)"
rm -rf plt0*
mpirun -np $NP "$EXE" inputs_couple > run_couple.log 2>&1 || { echo "run failed (see run_couple.log)"; exit 1; }
python3 check_sensible.py budget run_couple.log plt_diag_00040 plt00040 0.125 $CP || status=1
echo "== inflow ($NP ranks)"
rm -rf plt0*
mpirun -np $NP "$EXE" inputs_inflow > run_inflow.log 2>&1 || { echo "run failed (see run_inflow.log)"; exit 1; }
python3 - <<'PY' || status=1
import check_sensible as c, numpy as np
d = c.read_plotfile("plt00080"); th = d["theta"]; k = 4   # z = 22.5 m
up = th[40:50, 60:68, k].mean(); wake = th[70:90, 60:68, k].mean()
ok = wake > up + 0.01
print(f"  inflow/outflow: theta upstream {up:.3f} K, in the wake {wake:.3f} K -> {'PASS' if ok else 'FAIL'}")
raise SystemExit(0 if ok else 1)
PY
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
