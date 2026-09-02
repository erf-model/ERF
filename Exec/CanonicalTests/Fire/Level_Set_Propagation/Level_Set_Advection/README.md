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
| `inputs_fire_levelset_baseline` | Passing. Advection with signed-distance reinitialization. |
| `inputs_fire_levelset_no_reinit` | Passing. Control: advection only. |

## The level-set path uses a true signed distance

`phi` is in **metres** for this path — negative inside the burned region,
positive outside, `|grad phi| = 1`. The FARSITE path keeps its own normalized
`[-1, 1]` indicator convention, and `initialize_ignition` takes a flag to choose
between them.

This matters. A level-set method's advection, its Godunov Hamiltonian and the
Russo-Smereka subcell term are all derived for a signed distance in metres.
Running the solver on a `[-1, 1]`-clamped field flattens everything outside the
band, so the front eventually advances into ground carrying no gradient
information. That produced discontinuous jumps in burned area — 44 to 528 cells
in a single step — and widening the band made it worse rather than better, which
is what ruled out band width as the cause and indicted the normalization itself.

## Expected Results

Burned-cell count at t = 150 / 300 / 450 / 600 s, measured on 1 rank:

| | 150 s | 300 s | 450 s | 600 s |
|---|---|---|---|---|
| analytic, `r = r0 + R*t` | 34 | 41 | 48 | 56 |
| `inputs_fire_levelset_baseline` | 36 | 52 | 60 | 64 |
| `inputs_fire_levelset_no_reinit` | 32 | 44 | 52 | 52 |
| FARSITE reference | 32 | 32 | 32 | 68 |

Both cases track the analytic front. `phi` behaves as a signed distance:
`phi_min` about -38 m near the centre of the burn, `phi_max` about the
far-corner distance, unclamped.

The table above is measured at `OMP_NUM_THREADS=1`.

### Determinism

Repeat runs agree exactly, and both cases are clean under
`amrex.init_snan=1 amrex.fpe_trap_invalid=1`.

Neither was true before the fire `Geometry` was given the atmospheric
periodicity. `create_fire_grid()` hard-coded `{false, false, false}`, so every
`FillBoundary` on the fire grid was a no-op: the fire grid is a single box
spanning the domain, which makes all of its ghost cells domain-boundary ghosts.
The WENO5-Z stencil reaches three cells past the box edge and was reading
whatever the allocator supplied, so the burned area varied run to run — 64, 82,
106, 118 and 126 cells at 600 s across five identical runs.

The FARSITE path was unaffected, which matches operational experience with
periodic atmospheres. It rebuilds `phi` from `fire_arrival_time` every step so
nothing accumulates, and its neighbour access is guarded with in-box bounds
(`i+1 <= hi.x`) rather than reaching into ghosts. It exits cleanly under the
signaling-NaN trap that made the level-set path abort.

## Defects fixed in the reinitializer

| Defect | Symptom before the fix |
|---|---|
| `phi` clamped to `[-1, 1]` rather than a signed distance in metres | Front advanced into a flat field; burned area jumped 44 to 528 cells in one step |
| Update targeted `\|grad phi\| = 1` while `phi` was normalized | `phi` diverged to about `1e7` in 20 steps; the jump equalled `n_iters * dtau` exactly |
| Gradient took the larger-magnitude one-sided difference, not the Godunov upwind | Burned cells flooded from 32 to 20786, over half the domain |
| No subcell fix, so the iteration moved the zero level set | Each pass eroded the front: 32 to 24 cells, and the fire went out entirely by 100 s |
| Default `dtau = 0.5*dx` sat exactly on the `dtau <= dx/2` limit | Unstable at 10 or more iterations |

The Russo-Smereka (2000) subcell update now fixes the interface from `phi_0` for
cells whose original neighbourhood straddles it, so a pass no longer erodes the
front: 32 cells before and after, at every iteration count tested from 5 to 80.

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
