# SEB/Phase8_WallFunction

Phase 8, first part: the wall function beyond neutral. Two switches, both
off by default, so every earlier result is unchanged unless a deck asks
for them.

```
./run_wallfunction.sh /path/to/erf_exec        # NP=4 by default, a few minutes
```

## The problem

The neutral log law on the tangential wind cannot shed heat from a hot
face in calm air: with the wind at its 1e-3 m/s floor the friction
velocity vanishes and so does the flux. In the day canonical
(`Phase7_IsolatedBuilding`) the roof reaches 340 K largely for that
reason.

## The switches

- `erf.ibseb.convective_velocity = deardorff`: the wind the wall function
  sees becomes `sqrt(U_tan^2 + (beta w*)^2)` (Beljaars' gustiness form)
  with `w* = (g/theta H/(rho c_p) depth)^(1/3)` from the previous step's
  flux out of the face (zero when the flux is into it). The depth is the
  mixed layer above a roof, `z_i - z_face` floored at the building
  height, and the building height for a wall. `z_i` follows the day:
  `erf.ibseb.z_i_mode = bulk_ri` diagnoses it every step on the
  horizontal-mean profile (first cell centre where the bulk Richardson
  number exceeds `ri_crit`; the domain top on a neutral profile), `pblh`
  reads the surface layer's boundary-layer height at the column, `fixed`
  takes `erf.ibseb.z_i`.
- `erf.ibseb.stability_correction = true`: Dyer's similarity functions on
  the roofs, iterated to convergence on the face's own Obukhov length
  `L = u*^2 theta / (kappa g theta*)`, under-relaxed as the surface
  layer's iteration is in erf-model #3486, and seeded from the ground
  surface layer's 2D field at the face's column. Walls stay on the log
  law: the functions assume a horizontal surface, and on a wall the
  convective scale carries free convection.

## The scenario

The phase 7 cube (an exact box) in calm air under a fixed sun 10 degrees
off the zenith at 1000 W/m2, with a light 5 cm cladding so the roof heats
within minutes; 600 steps of 0.5 s, every step dumped. The ground surface
layer is neutral (300 K on 300 K air), so the seed is neutral and the
roofs' L is entirely their own.

## The decks and what is checked

1. `inputs_neutral` against `inputs_deardorff` (fixed 1000 m depth): the
   roof at 341 K sheds under 1 W/m2 with the neutral law and hundreds
   with the convective scale.
2. `inputs_deardorff`: on every face of the last dump the depth is the
   mixed layer above the roofs and the building height on the walls, w*
   is the Deardorff scale of the previous dump's flux, and u* and H
   follow the log law on the effective wind, all to 1e-9.
3. `inputs_stability` (functions and the scale): the roofs' Obukhov
   length is negative and equals `u*^2 theta / (kappa g theta*)` with
   theta* from the previous step's skin (what the wall function saw; the
   stored H is the balance's flux at the new skin), u* and H follow the
   corrected log law with Dyer's psi_m and psi_h at delta/L to 1e-7, and
   the roof flux exceeds the run without the functions.
4. `inputs_bulkri` (capped sounding, inversion at 100 m): the diagnosed
   depth is the first cell centre above the inversion, 95 m, and the
   roofs' depth in w* is that minus the 40 m roof.

## Reference output (4 ranks)

```
== neutral (4 ranks, 600 steps)
T_skin_min=304.6169719 T_skin_max=341.5596623 SW_abs_mean=340.6282257 SW_abs_max=976.3269777 shadow_frac=0 LW_net_mean=-137.1782516 H_mean=0.06753376039 G_mean=203.3824403 Q_ext_mean=0 resid_max=2.242131814e-09 w_star_max=0 H_total_W=540.270083109
== deardorff (4 ranks, 600 steps)
T_skin_min=304.4643253 T_skin_max=331.6713707 SW_abs_mean=340.6282257 SW_abs_max=976.3269777 shadow_frac=0 LW_net_mean=-119.36381 H_mean=59.22946039 G_mean=162.0349554 Q_ext_mean=0 resid_max=2.171532287e-09 w_star_max=1.941227363 H_total_W=473835.683097
== stability (4 ranks, 600 steps)
T_skin_min=304.4754743 T_skin_max=324.7544115 SW_abs_mean=340.6282257 SW_abs_max=976.3269777 shadow_frac=0 LW_net_mean=-108.9114065 H_mean=92.77441933 G_mean=138.9423999 Q_ext_mean=0 resid_max=2.127933385e-09 w_star_max=2.281911818 H_total_W=742195.354618
== bulkri (4 ranks, 600 steps)
T_skin_min=304.456798 T_skin_max=338.3709009 SW_abs_mean=340.6282257 SW_abs_max=976.3269777 shadow_frac=0 LW_net_mean=-130.1564219 H_mean=23.9846172 G_mean=186.4871866 Q_ext_mean=0 resid_max=2.223373485e-09 w_star_max=0.5254213395 H_total_W=191876.937604
  roof in calm air: neutral law sheds little, the convective scale ten times more: PASS (roof H neutral 0.18 W/m2 at T_skin 341.6 K, with w* 267.8 W/m2 at 331.4 K; U_tan 1.5e-03 m/s)
neutral: PASS
  depth in w*: mixed layer above the roofs, building height on the walls: PASS (roofs 960 m, walls 40 m)
  w* is the Deardorff scale of the previous step's flux: PASS (max diff 6.0e-12 m/s, w* up to 1.94 m/s)
  u* and H from the log law on sqrt(U_tan^2 + (beta w*)^2): PASS (rel diff u* 3.6e-12, H 2.5e-11)
deardorff: PASS
  roof Obukhov length negative (unstable): PASS (L in [-2.7, -2.6] m)
  L = u*^2 theta / (kappa g theta*) with the face's own u* and theta* (previous skin): PASS (rel diff 3.7e-07)
  u* and H from the log law with Dyer's psi_m, psi_h at delta/L on the roofs: PASS (rel diff u* 7.1e-08, H 1.5e-07)
  roof flux above the run without the functions: PASS (with 436.1 W/m2, without 267.8 W/m2)
stability: PASS
  bulk Richardson depth on the capped sounding: PASS (diagnosed 95.0 m (expected 95), range over the run 95.0-95.0 m)
  roofs' depth in w* is the mixed layer above the roof: PASS (55.0 m above the 40 m roof)
bulkri: PASS
ALL PASS
```
