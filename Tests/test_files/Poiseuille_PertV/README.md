# Poiseuille_PertV

## Purpose

Checks, for erf-model/ERF#4618, that the y-velocity perturbation of the
Couette/Poiseuille setup (`ERF_InitCustomPertVels_CouettePoiseuille.H`) follows
`prob.pert_delta_v` and only that input.

## Files

| file | role |
| --- | --- |
| `Poiseuille_PertV.i` | the channel setup (`prob.prob_type = 21`) on 16 x 8 x 16 cells, initial state only, `prob.pert_delta_u = 0`, `prob.pert_delta_v = 2` |
| `../../check_v_perturbation.py` | the check (Python standard library only) |

## How to run

```bash
mpiexec -n 2 erf_exec Poiseuille_PertV.i
python3 ../../check_v_perturbation.py plt00000 --expect nonzero
mpiexec -n 2 erf_exec Poiseuille_PertV.i prob.pert_delta_u=2.0 prob.pert_delta_v=0.0
python3 ../../check_v_perturbation.py plt00000 --expect zero
```

## Expected results

The setup has no mean y-velocity, so `y_velocity` in the initial plotfile is the
perturbation alone.

| run | before the fix | this code |
| --- | --- | --- |
| `pert_delta_v = 2` only | max abs(v) = 0 | 1.84 |
| `pert_delta_u = 2` only | 1.84 | 0 |

Before the fix `pert_delta_v` was read from the key `prob.pert_delta_u`, and the
y-velocity block was switched on by `pert_delta_u`.

## CTest

`Poiseuille_PertV_Only` and `Poiseuille_PertU_Only` (label `regression`; need Python 3).
