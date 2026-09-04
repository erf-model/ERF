# ImmersedForcingTest/PartialCells

Regression test for the immersed forcing on the partial cells of a
height-map building, and for the input `erf.if_snap_partial_cells`.

```
./run_partial_cells.sh /path/to/erf_exec               # NP=4, about 10 minutes
./run_partial_cells.sh /path/to/erf_exec --reproduce   # also runs the original selection, which traps
```

## The problem

The terrain reader interpolates a nodal height map, so every building edge
becomes a one-cell ramp and its corners are sliver cells with 1 to 20
percent solid. The forcing has two regimes on a partial cell: a drag
proportional to the solid fraction, which fades to zero with it, and the
wall law (with `erf.if_use_most`), which the original selection applies
to any cell above 0.005 solid whose normal neighbour is exactly zero, at a
rate that does not depend on the fraction. The top cell of a sliver
column is therefore forced like a full wall on top of an almost free
cell. In a neutral 3 m/s run on a 40 m cube at 10 m and 0.5 s, a vertical
two-cell checkerboard in theta grows there over hours, reaching 10 K with
5 m/s vertical velocity spikes, until the floating-point trap fires in
the Exner function after 1.8 to 2.5 hours. It does so with the surface
energy balance off, with a no-slip ground, without the substep forcing
(negative density at 3 minutes) and with `eb2.small_volfrac = 0.25` (trap
at step 14); halving the step only delays it.

## The switch

`erf.if_snap_partial_cells = true` makes the six forcing functions agree
on where the wall is: the wall law and the surface temperature, flux and
Obukhov conditions sit on cells at least half solid whose normal
neighbour is less than half solid, and cells below half solid keep only
the fraction-weighted drag. The rule is the one the surface energy balance
uses to find its faces, so the two coincide. Under the switch a face at
a ramp or a corner can pass several wall tests at once (a half-blanked
roof-edge face sees its side neighbours at a quarter, which now count as
fluid), so the stacked relaxations are averaged: the total rate on a face
never exceeds that of a single wall, which the explicit forcing in the
substeps needs (without this the density goes negative in three steps).
The default is false, the original behaviour, for backward
compatibility; buildings from height maps should set it. With the switch
off every result is bit-identical to before.

## What is checked

The run with the switch completes 2.64 hours (the original traps), the
fluid temperature stays within 0.5 K of neutral, the vertical velocity
below 2 m/s, and the cube still throws a wake.

## Reference output (4 ranks)

```
== height-map cube, wall law snapped to half cells (4 ranks, 19000 steps)
  fluid temperature within 0.5 K of neutral: PASS (max |theta - 300| 0.000 K at t = 9500 s)
  vertical velocity below 2 m/s: PASS (max |w| 0.49 m/s)
  wake behind the cube: PASS (min u in the wake -0.38 m/s (inflow 3))
partial cells: PASS
ALL PASS
```
