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

### Open defect: the reinitialized result depends on thread count

| Threads | 150 s | 300 s | 450 s | 600 s |
|---|---|---|---|---|
| 1 | 36 | 52 | 60 | 64 |
| 2 | 36 | 52 | 62 | 98 |
| 4 | 36 | 52 | 86 | 118 |

Repeatable within a thread count, different across. `inputs_fire_levelset_no_reinit`
is thread-independent (32/44/52/52 at both 1 and 4 threads), which localises the
defect to the reinitialization rather than the advection. Compare against the
single-thread column until it is resolved.

A separate non-determinism — differing run to run at a fixed thread count — was
traced to the Jacobi scratch fab being left uninitialized, so the result depended
on whatever the allocator handed over; that one is fixed.

A CPU threading defect of this kind usually has a GPU counterpart, so this is
worth resolving before the level-set path is trusted on device.

**Reinitialization must be converged.** It propagates information roughly
`n_iters * dtau` metres per call, so the default 10 iterations only reaches
about 25 m here and leaves the field short of a signed distance further out.
That shows up as the front over-spreading — 122 cells at 600 s against 64 with
100 iterations. Scale `reinit_iters` to the size of the burned region, not to
the cell size.

Two further properties are worth checking because they are not visible in a
single-rank run:

- **No box-boundary artefacts in `phi`.** The RK3 stages fill ghost cells before
  every WENO5-Z reconstruction, whose stencil reaches three cells past a box
  edge. Run on several ranks and look for structure aligned with box boundaries.
- **Bitwise reproducibility.** Both the advection and the reinitialization read
  neighbours while writing the centre cell, so both are done Jacobi-style into
  scratch. Repeat runs should agree exactly, on CPU and GPU alike.

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
