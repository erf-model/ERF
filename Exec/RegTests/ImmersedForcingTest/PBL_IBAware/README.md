# ImmersedForcingTest/PBL_IBAware

Regression test for `erf.pbl_ib_aware`: the MRF and YSUNew boundary-layer
schemes measuring height from the top of an immersed building instead of
the terrain surface.

```
./run_pbl_ib.sh /path/to/erf_exec        # NP=4, a few minutes
```

## The problem

The PBL schemes take the wall distance from the terrain surface (or the
domain bottom on a flat grid), so over an immersed building the mixing
length `kappa z (1 - z/h)^2`, the bulk Richardson heights and the
boundary-layer depth are all off by the building height, the diffusivity
profile continues inside the solid, and the surface scales u*, theta*, L
of the column come from the ground surface layer evaluated on cells
inside the building.

## The switch

With `erf.pbl_ib_aware = true` (MRF and YSUNew only; the other schemes are
untouched) each column's surface is its first fluid cell above the
immersed solid, found from the blanking every step: every height in the
Richardson diagnosis and the K profile is measured from it, the
boundary-layer depth with it, the diffusivities are zero inside the
solid, and the surface scales of a column with solid cells are a neutral
log law at its top (u* from the wind there and `erf.pbl_ib_z0`, theta* = 0,
L large). Columns without solid cells are untouched, so a domain without
immersed cells gives bit-identical results. Not supported with
terrain-fitted coordinates. Default false.

## What is checked

On a 40 m box cube in a 3 m/s neutral flow after 100 s, for each scheme
with the switch on:

1. The vertical eddy viscosity is finite everywhere and zero in every
   solid cell.
2. The roof column's profile against height above the roof follows the
   open-ground column's against height above the ground: the same shape,
   and the same peak height (15 m for MRF; YSUNew is flat in this neutral
   case and only the shape is compared). The wall distance restarts at
   the roof.
3. The first cell above the roof is within a factor five of the first
   cell above the ground.
4. Without a building the switch changes nothing (bit-identical).

The runs with the switch off are reported, not asserted: on this case MRF
completes with the viscosity NaN over most of the domain (the ground
surface layer evaluates u* and L on cells inside the building and the
NaN spreads), and YSUNew drives the density negative at its second step.
That is the defect the switch removes; without immersed cells both
schemes are untouched. MRF also has an invalid floating-point operation
at its neutral start with or without a building, so its decks run with
the trap off.

## Reference output (4 ranks)

```
== mrf_off (4 ranks, 200 steps; may fail, see the README)
== ysunew_off (4 ranks, 200 steps; may fail, see the README)
   run with the switch off failed, as documented
== mrf_on (4 ranks, 200 steps)
== ysunew_on (4 ranks, 200 steps)
== flat_off (4 ranks, 200 steps)
== flat_on (4 ranks, 200 steps)
  switch on: finite everywhere, zero in every solid cell: PASS (0 non-finite, max K inside the cube 0.0e+00)
  switch on: the roof column's profile against height above the roof follows the ground's against height above the ground: PASS (shape difference 0.23; peak 15 m above the roof, 15 m above the ground)
  switch on: first cell above the roof within a factor five of the first cell above the ground: PASS (roof 0.625, ground 0.138)
  switch off: run completed with 16320 non-finite viscosity cells (unusable)
scheme: PASS
  switch on: finite everywhere, zero in every solid cell: PASS (0 non-finite, max K inside the cube 0.0e+00)
  switch on: the roof column's profile against height above the roof follows the ground's against height above the ground: PASS (shape difference 0.00; flat profile)
  switch on: first cell above the roof within a factor five of the first cell above the ground: PASS (roof 0.127, ground 0.128)
  switch off: run failed: Rho is negative at 16 16 5 -4.618879e+00
scheme: PASS
  no building: switch on and off bit-identical: PASS (max |dK| 0.0e+00)
flat: PASS
ALL PASS
```
