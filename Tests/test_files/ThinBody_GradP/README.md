# ThinBody_GradP

## Purpose

Checks, for erf-model/ERF#4608, that the thin-body projection (`project_velocity_tb`) keeps the pressure
gradient `gradp` consistent with the pressure increment `pp_inc` when it makes more
than one corrector pass (`erf.ncorr > 1`) and when dt is not 1.

## Files

| file | role |
| --- | --- |
| `ThinBody_GradP.i` | anelastic flow in a closed 2 x 2 x 1 box, 16 x 16 x 8 cells, 2 x 2 x 1 boxes, a thin wall of 16 zero-x-flux faces, `erf.ncorr = 2`, dt = 0.025, 2 steps |
| `../../check_thin_body_gradp.py` | the check (Python standard library only) |

## How to run

```bash
mpiexec -n 2 erf_exec ThinBody_GradP.i
python3 ../../check_thin_body_gradp.py plt00002
```

## Expected results

The cell average of the two face gradients around a cell is the centred difference
of the cell values, so in every interior cell `dpdx` equals
`(pert_pres(i+1) - pert_pres(i-1)) / (2 dx)`, and the same in y and z. The check
passes when the largest mismatch is below 1e-8 of the largest gradient; the fixed
code gives about 3e-16. It stops if the gradient is zero everywhere.

The check compares two stored fields with each other, so a wrong factor common to both
would pass; the 1/dt in `pp_inc` is pinned by `gradp`, which already carried it
(the same convention as `project_momenta`). Double precision only.

| build | relative mismatch |
| --- | --- |
| before the fix (`gradp` from the last pass only, `pp_inc` without 1/dt) | 1.0 |
| `pp_inc` fixed, `gradp` still from the last pass only | 12.9 |
| this code | 3e-16 |

## CTest

`ctest -R ThinBody_GradP` (label `regression`; needs Python 3).
