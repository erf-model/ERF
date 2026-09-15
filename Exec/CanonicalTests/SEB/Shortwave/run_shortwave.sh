#!/bin/bash
# Regression test: shortwave with ray-cast shadowing.
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
hdr = open("faces_zen60_np.csv").readline().strip().split(",")
a = load("faces_zen60_np.csv"); b = load("faces_zen60.rank0.csv")
# Geometry, view fractions, shadow, shortwave and materials must agree exactly;
# the columns read from the atmosphere (temperatures, wind, fluxes, skin) may
# differ by the round-off of the decomposition, about 1e-12 relative after
# two steps, which the dump's tenth digit shows.
exact = {"i", "j", "k", "dir", "side", "bid", "x_m", "y_m", "z_m", "area_m2", "f_sky", "f_ground", "f_bldg", "shadow",
         "SW_direct_in", "SW_diffuse_in", "SW_abs", "mat", "albedo", "emissivity", "k_therm", "rho_cp", "thickness", "h_bld"}
bad = [] if a.shape == b.shape else ["shape"]
for n, h in enumerate(hdr):
    if a.shape == b.shape and not np.allclose(a[:, n], b[:, n], rtol=0 if h in exact else 1e-9, atol=0 if h in exact else 1e-9): bad.append(h)
ok = not bad
print(f"  rank independence of the face dump: {'PASS' if ok else 'FAIL'} ({a.shape[0]} vs {b.shape[0]} faces; geometry, view, shadow, shortwave exact, atmosphere columns to 1e-9{', differing: ' + ' '.join(bad) if bad else ''})")
raise SystemExit(0 if ok else 1)
PY
[ $status -eq 0 ] && echo "ALL PASS" || echo "SOME CHECKS FAILED"
exit $status
