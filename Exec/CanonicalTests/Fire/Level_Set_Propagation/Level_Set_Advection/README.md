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

## Expected Results
Both cases should complete with a monotonically growing burned area and `phi`
bounded in `[-1, 1]`. The two differ only in how hard the reinitializer is
driven, so their burned areas should stay close; a large divergence points at
the reinitialization rather than the advection.

Two properties are worth checking explicitly because they are not visible in a
single-rank run:

- **No box-boundary artefacts in `phi`.** The RK3 stages fill ghost cells before
  every WENO5-Z reconstruction, whose stencil reaches three cells past a box
  edge. Run on several ranks and look for structure aligned with box boundaries.
- **Bitwise reproducibility.** The reinitialization update reads neighbours while
  writing the centre cell, so it is done Jacobi-style into scratch. Repeat runs
  should agree exactly, on CPU and GPU alike.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.fire.propagation_method` | `"levelset"` | Selects the PDE solver over the FARSITE Lagrangian default. |
| `erf.fire.levelset.cfl` | `0.4` / `0.25` | Subcycle CFL number. Must be `> 0`; a non-positive value is rejected at startup. |
| `erf.fire.levelset.eps_visc` | `0.4` / `0.2` | Artificial viscosity coefficient on the Laplacian term. |
| `erf.fire.levelset.reinit_every` | `5` / `1` | Reinitialize every N subcycles. Must be `>= 1`; it is a modulus divisor. |
| `erf.fire.levelset.reinit_iters` | `10` / `20` | Sussman pseudo-time iterations per reinitialization. |
| `erf.fire.levelset.reinit_dtau` | `-1.0` / `2.0` | Pseudo-time step; `<= 0` selects `0.5*min(dx,dy)`. |

## References
- Osher & Sethian 1988, Fronts propagating with curvature-dependent speed.
- Sussman, Smereka & Osher 1994, A level set approach for computing solutions to incompressible two-phase flow.
- Borges et al. 2008, An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws (WENO-Z).
