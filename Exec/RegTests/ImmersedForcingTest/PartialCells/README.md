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

`erf.if_snap_partial_cells = true` makes the six forcing functions read the
blanking snapped to solid or fluid at half: a height-map building becomes
the same staircase of whole cells an exact box is, the wall law sits on
the faces between a solid and a fluid cell, the drag inside the solid, and
the sliver cells carry nothing. The 24 h day of `SEB/Phase7_IsolatedBuilding`
on an exact box is the evidence that this configuration is stable. (A
first version that kept the fractions and only moved the selection
threshold to half left faces between a solid core and a half-solid rim
relaxed toward a wall-law target while the core's other faces were damped,
and a thin slab drove its density negative within two minutes.) The default
is false, the raw fractions, for backward compatibility; buildings from
height maps should set it. With the switch off every result is
bit-identical to before, and on an exact box the switch changes nothing.

## What is checked

The run with the switch completes 2.64 hours (the original traps), the
fluid temperature stays within 0.5 K of neutral, the vertical velocity
below 2 m/s, and the cube still throws a wake.

## Reference output (4 ranks)

```
== height-map cube, wall law snapped to half cells (4 ranks, 19000 steps)
  fluid temperature within 0.5 K of neutral: PASS (max |theta - 300| 0.000 K at t = 9500 s)
  vertical velocity below 2 m/s: PASS (max |w| 0.60 m/s)
  wake behind the cube: PASS (min u in the wake -0.32 m/s (inflow 3))
partial cells: PASS
ALL PASS
```
