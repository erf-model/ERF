#!/bin/bash
# Regression test: view fractions and longwave.
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
hdr = open("faces_fixed_np.csv").readline().strip().split(",")
a = load("faces_fixed_np.csv"); b = load("faces_fixed.rank0.csv")
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
