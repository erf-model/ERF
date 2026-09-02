# Level_Set_Advection

## Purpose
Exercises the level-set fire-front solver: WENO5-Z gradient reconstruction, a
three-stage SSP-RK3 step subcycled on its own CFL condition, and Sussman
signed-distance reinitialization.

## Physics / Model Features Exercised
- Level-set advection of the fire front (`propagation_method = "levelset"`)
- WENO5-Z reconstruction with artificial viscosity on the Laplacian
- CFL-based subcycling within one atmospheric step
- Sussman reinitialization of the signed-distance property

## Status

| Case | State |
|---|---|
| `inputs_fire_levelset_baseline` | Passing. Advection only; reinitialization disabled. |
| `inputs_fire_levelset_reinit_erosion` | **Known limitation.** Reproducer for front erosion under reinitialization. |

## Expected Results

**Baseline.** Completes with a monotonically growing burned area and `phi` in
`[-1, 1]`. Measured at the 100 s stop time on 1 rank: `phi` in `[-0.84, 1.0]`,
34 burned cells. For reference the FARSITE path holds `phi` in `[-1, 0]`.

Two further properties are worth checking explicitly because they are not
visible in a single-rank run:

- **No box-boundary artefacts in `phi`.** The RK3 stages fill ghost cells before
  every WENO5-Z reconstruction, whose stencil reaches three cells past a box
  edge. Run on several ranks and look for structure aligned with box boundaries.
- **Bitwise reproducibility.** Both the advection and the reinitialization update
  read neighbours while writing the centre cell, so both are done Jacobi-style
  into scratch. Repeat runs should agree exactly, on CPU and GPU alike.

## Reinitialization: fixed defects and the one that remains

Three defects were repaired in the reinitializer:

| Defect | Symptom before the fix |
|---|---|
| Update targeted `\|grad phi\| = 1`, i.e. a signed distance in **metres**, but ERF-Hazard's `phi` is normalized to `[-1, 1]` | `phi` diverged to about `1e7` within 20 steps; the jump equalled `n_iters * dtau` exactly |
| Gradient took the larger-magnitude one-sided difference instead of the Godunov upwind | Burned cells flooded from 32 to 20786, over half the domain |
| Default `dtau = 0.5*dx` sat exactly on the `dtau <= dx/2` stability limit | Unstable at 10 or more iterations |

The update now targets `|grad phi| = 1/L` for an explicit band half-width `L`
(`erf.fire.levelset.reinit_band_m`, default `3*min(dx,dy)`), selects the
one-sided difference by the sign of `phi_0`, and defaults `dtau` to `0.25*dx`.

**What remains.** The Sussman iteration still does not preserve the zero level
set, so a reinitialization pass erodes the front:

| Iterations | Burned cells across one pass |
|---|---|
| 5 | 32 -> 24 |
| 10 | 32 -> 24 |
| 20 | 32 -> 24 |
| 40 | 32 -> 16 |

Interface displacement is a known property of the basic scheme and is normally
handled with a subcell fix (Russo & Smereka 2000), which constrains cells
adjacent to the interface using the pre-reinitialization values. That is not
implemented here. The loss is permanent: unlike the FARSITE path, which rebuilds
`phi` from `fire_arrival_time` every step, the level-set path has no
reconstruction, so an eroded front never recovers.

**Consequence.** The baseline runs with reinitialization disabled and stops at
100 s. Without reinitialization `phi`'s magnitude drifts — measured `phi_min` is
`-0.84` at 100 s, `-1.57` at 300 s, `-3.40` at 600 s — while its sign, which is
all the rest of the module reads, stays correct and the front keeps advancing.
Until the subcell fix lands, treat the level-set path as experimental and prefer
the FARSITE default for production runs.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.fire.propagation_method` | `"levelset"` | Selects the PDE solver over the FARSITE Lagrangian default. |
| `erf.fire.levelset.cfl` | `0.4` / `0.25` | Subcycle CFL number. Must be `> 0`; a non-positive value is rejected at startup. |
| `erf.fire.levelset.eps_visc` | `0.4` / `0.2` | Artificial viscosity coefficient on the Laplacian term. |
| `erf.fire.levelset.reinit_every` | `1000000` / `1` | Reinitialize every N subcycles. Must be `>= 1`; it is a modulus divisor. Set high in the baseline to disable it. |
| `erf.fire.levelset.reinit_iters` | `10` / `20` | Sussman pseudo-time iterations per reinitialization. |
| `erf.fire.levelset.reinit_dtau` | `-1.0` / `2.0` | Pseudo-time step; `<= 0` selects `0.5*min(dx,dy)`. |

## References
- Osher & Sethian 1988, Fronts propagating with curvature-dependent speed.
- Sussman, Smereka & Osher 1994, A level set approach for computing solutions to incompressible two-phase flow.
- Borges et al. 2008, An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws (WENO-Z).
