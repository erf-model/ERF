# Canonical RANS cases

Regression cases for the one-equation k RANS closure of Axell & Liungman
(2001), `erf.rans_type = kEqn`, laid out like `../Canonical_LES`: one
directory per case with the input deck, the sounding, a README stating the
targets and tolerances, and a Python check script. The theory and the
inputs are described in `Docs/sphinx_doc/theory/RANS.rst`; the development
record is `PLAN.md` and the numbers per phase are in `RESULTS.md`.

| case | mesh | physics run | what it checks |
| --- | --- | --- | --- |
| `Neutral_ABL_Flat` | flat, 8 x 8 x 64 | 12 h | log law, wall k = u*^2/Cmu0^2, length scale, dissipation |
| `Stable_ABL_Flat` | flat, 8 x 8 x 100 | 9 h (GABLS1) | u*, low-level jet, depth, stable stratification |
| `Convective_ABL_Flat` | flat, 8 x 8 x 100 | 4 h | column heat budget, inversion height, mixed layer, wall k with buoyancy |
| `Neutral_Hill_2D` | fitted, 128 x 1 x 64 | 6 h | wall distance vs the exact ridge distance, crest speed-up, upstream log law |
| `Neutral_Hill_3D` | fitted, 64 x 64 x 20 | 4 h | wall distance vs the exact hill distance, crest speed-up, upstream log law |

## Rules

- Every deck ships a `check_<case>.py` that reads the plotfile
  (`erf_plotfile.py`, standard library only) and prints one row per check:
  measured value, target, tolerance, pass. Its exit code is the verdict.
  `--smoke` runs the structural checks that must hold after a few steps
  (the CTest entries, `ctest -L rans`); `--physics` adds the checks that
  need the converged run in the table above.
- Shared checks live in `rans_checks.py`; the terrain scripts use the full
  3D reader because ERF writes no planar averages on a fitted mesh.
- Closure changes are cross-checked against the Kynema `KLAxell` and
  `KransAxell` implementation and the paper; the comparisons are recorded
  in `PLAN.md`.

## Running a case

```bash
cd Neutral_ABL_Flat
mpirun -np 2 erf_exec inputs_neutral
python3 check_neutral.py --physics plt08640 surf_hist.dat
```

`surf_hist.dat` (u*, theta*, L) is written only when `erf.v > 0`.
