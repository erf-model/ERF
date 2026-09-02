# Level_Set_Propagation

## Purpose
This directory collects canonical cases for the level-set fire-front propagation
method, the PDE-based alternative to the default FARSITE Lagrangian path selected
with `erf.fire.propagation_method = "levelset"`.

## Subcases
| Subdirectory | Description | Character |
|---|---|---|
| `Level_Set_Advection` | WENO5-Z advection with SSP-RK3 time stepping and Sussman signed-distance reinitialization. | empirical / regression |

## Notes
- Level-set propagation is **isotropic**. The front advances at the scalar rate
  of spread through a Godunov norm and does not use wind direction, so a
  wind-driven fire stays near-circular rather than developing the ellipse the
  FARSITE path produces. Compare against `FARSITE_Propagation/`.
- These cases are the only coverage of `erf.fire.levelset.*`, so they should be
  run whenever the level-set solver changes.
