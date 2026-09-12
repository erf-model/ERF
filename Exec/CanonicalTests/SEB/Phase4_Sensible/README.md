# SEB/Phase4_Sensible

Phase 4 of the immersed-boundary surface energy balance: sensible heat from
every face through a wall function, added to the temperature equation of
the adjacent fluid cell. The skin temperature is prescribed; it does not
evolve until phase 6.

```
./run_sensible.sh /path/to/erf_exec        # NP=4 by default
python3 plot_sensible.py plt00080          # writes plots/*.png
```

## The scenario

A 40 m cube at the centre of a 640 m domain on a 5 m grid in a uniform
8 m/s wind at 300 K, its faces held at 320 K. No sun, no Coriolis, no
damping, an adiabatic no-slip ground, so the only heat entering the air is
the face flux. The immersed forcing carries the momentum wall model
(`erf.if_use_most`); its own surface-temperature inputs must be absent, the
balance owns the temperature condition at the buildings.

## The wall function

On each face the tangential wind of the fluid cell at half a cell from the
wall gives the friction velocity with the roughness `z0_wall`, and the
skin-to-air potential-temperature difference gives the temperature scale
with `z0h_wall`; H = rho c_p u* theta* is positive out of the face. Neutral
on every face; `erf.ibseb.stability_correction` applies the surface
layer's similarity functions on roofs. The flux enters the rho-theta
equation as H A / (c_p V Pi) with Pi the Exner function, so the fluid cell
warms by H A / (c_p V) per second in temperature.

## The decks

- `inputs_diag`: the flux is diagnosed, not applied (`couple_heat = false`).
- `inputs_couple`: the flux is applied; 40 steps in the closed domain.
- `inputs_inflow`: mass inflow at the west with a uniform profile
  (`xlo.type = Inflow` with a dirichlet file, as the Askervein canonical),
  pressure outflow at the east, periodic in y; 80 steps.

## What is checked

1. On every face of the diagnostic dump, u* and H equal the wall-function
   formulas recomputed from the dumped wind, density and temperatures; H is
   positive everywhere and largest on the windward wall.
2. The diagnostic dump on one rank equals the four-rank dump.
3. With the flux applied, the heat the faces added over the run (sum of
   H_total dt from the log) matches the extra internal energy of the air
   against the diagnostic run at the same step within 10 percent. The
   difference between the two runs isolates the heating, since the whole
   domain's energy drifts by far more than the face heat while the initial
   state adjusts; and it is the internal energy, c_v = c_p - R_d, because
   the domain is a rigid closed box, so the enthalpy rises by c_p/c_v times
   the heat added (the ratio 1.37 was measured before this was accounted
   for).
4. The inflow/outflow deck runs and the wake behind the cube is warmer than
   the upstream air.

## Reference

```
== diag (4 ranks)
T_skin_min=320 T_skin_max=320 SW_abs_mean=0 SW_abs_max=0 shadow_frac=0 LW_net_mean=-358.4161197 H_mean=463.9740461 H_total_W=4500548.24729
  u* = kappa U_tan / ln(delta/z0) on every face        PASS 
  H = rho c_p u* kappa (theta_skin - theta_air)/ln(delta/z0h) PASS 
  H > 0 on every face (walls warmer than the air)      PASS H 4.8-941.8 W/m2
  H largest on the windward face (largest U_tan)       PASS 
== diag (1 rank)
  rank independence of the face dump: PASS (388 vs 388 faces)
== couple (4 ranks)
  heat from the faces 24.804 MJ, extra internal energy of the coupled air 24.247 MJ, relative difference 2.2% -> PASS
== inflow (4 ranks)
  inflow/outflow: theta upstream 300.000 K, in the wake 300.050 K -> PASS
ALL PASS

```
