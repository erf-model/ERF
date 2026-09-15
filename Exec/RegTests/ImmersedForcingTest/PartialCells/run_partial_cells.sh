#!/bin/bash
# Partial-cell regtest of the immersed forcing.
#
#   ./run_partial_cells.sh /path/to/erf_exec              # NP=4, about 10 minutes
#   ./run_partial_cells.sh /path/to/erf_exec --reproduce  # also run with the switch off (must trap)
set -u
EXE=${1:?usage: run_partial_cells.sh /path/to/erf_exec [--reproduce]}
NP=${NP:-4}
rm -rf plt*
echo "== height-map cube, wall law snapped to half cells ($NP ranks, 19000 steps)"
mpirun -np $NP "$EXE" inputs > run_snap.log 2>&1 || { echo "run failed (see run_snap.log)"; exit 1; }
python3 check_partial_cells.py plt19000 || { echo "SOME CHECKS FAILED"; exit 1; }
if [ "${2:-}" = "--reproduce" ]; then
    echo "== the same deck with the original selection (expected to trap within 2.5 h)"
    rm -rf plt*
    mpirun -np $NP "$EXE" inputs erf.if_snap_partial_cells=false > run_original.log 2>&1
    grep -c "SIGILL\|SIGFPE\|Bad rho" run_original.log > /dev/null && echo "  original selection trapped, as documented" || echo "  original selection did not trap this time (see run_original.log)"
fi
echo "== the wall law is live under the snap: z0 must change the flow (400 steps, $NP ranks)"
for z in 0.01 0.1; do
    rm -rf plt_z0_${z}_*
    mpirun -np $NP "$EXE" inputs max_step=400 erf.plot_int_1=400 erf.plot_file_1=plt_z0_${z}_ erf.if_z0=$z > run_z0_$z.log 2>&1 || { echo "run failed (see run_z0_$z.log)"; exit 1; }
done
python3 - <<'PY' || { echo "SOME CHECKS FAILED"; exit 1; }
import yt, numpy as np
yt.set_log_level(40)
a = yt.load("plt_z0_0.01_00400").all_data(); b = yt.load("plt_z0_0.1_00400").all_data()
d = np.abs(a["boxlib", "x_velocity"].value - b["boxlib", "x_velocity"].value).max()
ok = d > 1.0e-3
print(f"  z0 = 0.01 vs 0.1 after 400 steps: max |du| {d:.3e} m/s -> {'PASS' if ok else 'FAIL'} (a no-slip staircase gives 0)")
raise SystemExit(0 if ok else 1)
PY
echo "ALL PASS"
