# MassFlux_TrapInvalid

## Purpose

Checks the constant mass-flux control in `make_mom_sources` (`erf.const_massflux_u`,
`erf.const_massflux_v`, `erf.const_massflux_tau`, `erf.const_massflux_layer_lo/hi`)
with `amrex.fpe_trap_invalid = 1`.

Without a mass flux the cross-sectional areas of the control are zero. The control
divides by them only when they are positive, but an optimised build may divide first
and test the area afterwards. The 0/0 then trips the invalid-operation trap, although
the result is not used. Before the fix, Apple clang Release and RelWithDebInfo builds
stopped with SIGILL at step 1 (the code came in with erf-model/ERF#4749). Builds that
keep the test before the divide (for example Debug builds) pass either way.

## Files

| file | role |
| --- | --- |
| `MassFlux_TrapInvalid.i` | anelastic flow in a 2 x 1.5 x 0.8 m box, periodic in x and y, slip walls in z, 8 x 8 x 8 cells, dt = 0.025 s, `amrex.fpe_trap_invalid = 1` |
| `../../check_massflux_relaxation.py` | the relaxation check (Python standard library only) |

The sides differ, so the x and y cross-sections differ, and a control that used the
wrong one would fail the relaxation check.

## How to run

```bash
mpiexec -n 2 erf_exec MassFlux_TrapInvalid.i
mpiexec -n 2 erf_exec MassFlux_TrapInvalid.i max_step=3 erf.const_massflux_u=1.0 \
    erf.const_massflux_layer_lo=0.2 erf.const_massflux_layer_hi=0.6 > relax.log
python3 ../../check_massflux_relaxation.py relax.log --dt 0.025 --steps 3
```

## Expected results

**`MassFlux_Off_TrapInvalid`** runs one step without a mass flux and must finish.

**`MassFlux_X_Relaxation`** sets `erf.const_massflux_u = 1` m/s over the layer
0.2-0.6 m (cells 2-5, H = 0.4 m), with the default `erf.const_massflux_tau = 1` s. The
fluid starts at rest with density 1, so the target is T = 1 x 1 x 0.4 x 1.5 = 0.6 kg/s
and the integrated flux F starts at 0. Each RK2 step evaluates the control twice; with
a = dt / tau the flux after n steps is

    F_n = T (1 - r^n),   r = 1 - a + a^2 / 2

and the second evaluation of each step sees F_n + a (T - F_n). The check compares all six
printed values with these and requires the y flux to stay 0.

**`MassFlux_TauZero_Abort`** sets `erf.const_massflux_tau = 0` and must stop with
`erf.const_massflux_tau must be > 0`.

| build | `MassFlux_Off_TrapInvalid` | `MassFlux_X_Relaxation`, largest relative error | `MassFlux_TauZero_Abort` |
| --- | --- | --- | --- |
| development cdd97a697 (Apple clang Release) | SIGILL at step 1 | 1e-11 | fails (the run stops later in MLMG with a non-finite right-hand side) |
| this code | passes | 1e-11 | passes |
| this code, divide by the area left out | passes | fails | passes |
| this code, x divided by the y area | passes | fails (1.8e-2) | passes |

The 1e-11 comes from the ten digits of the printed values. The tolerance is 1e-8 in
double precision and 1e-4 in single precision.

## CTest

`ctest -R MassFlux_` (label `regression`; `MassFlux_X_Relaxation` needs Python 3).
