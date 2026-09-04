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
echo "ALL PASS"
