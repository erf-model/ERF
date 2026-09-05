# SEB/Phase6_Prognostic

Phase 6 of the immersed-boundary surface energy balance: the skin
temperature of every face is now solved from the balance every step,
consistently with the implicit slab, and the fluxes the air receives are
those of the closed balance. Phases 2-5 diagnosed each term around a fixed
skin; this phase closes the loop.

```
./run_prognostic.sh /path/to/erf_exec                       # NP=4 by default
python3 plot_prognostic.py faces_solar ibseb_solar.csv --plotfile plt10800 --dt 0.5
```

## The scenario

A 40 m cube at the centre of a 320 m periodic domain on a 10 m grid in a
3 m/s westerly at 300 K, its faces starting at 300 K over 20 cm of concrete
(k = 1 W/m/K, rho c_p = 1.6e6 J/m3/K, 8 layers) with a 293 K interior.
Gray sky, ground at 300 K, albedo 0.3, emissivity 0.9. The domain is small
so the 90 minute sunrise deck runs in a few minutes.

## The balance

Every step, after the shortwave, longwave and wall function of the step,
each face solves

    SW_abs + eps Q_ext + eps [LW_ext - (1 - f_bldg) sigma T^4] - C_H (T/Pi - theta_air) - LE - (a T - b) = 0

for T by Newton's method (`Source/ImmersedBoundarySEB/ERF_IBSEBBalance.H`,
lifted from the SLUCM facet solver). `a T - b` is the conduction the
implicit slab step will take for a skin temperature T, from two trial
slab steps, so the balance and the slab agree to rounding. `Q_ext` is an
incident external flux absorbed with the emissivity, the hook a fire's
radiation will use. The bounds `T_skin_min` / `T_skin_max` and the step cap
are inputs. After the solve the slab is advanced with the solution and
LW_net, H and G are rewritten at it.

## The decks

- `inputs_closure`: fixed sun in the south at 45 degrees, 200 steps of
  0.5 s, every step dumped (`dump_faces_tag_step`).
- `inputs_qext`: night, 3000 W/m2 incident on every face, a light thin
  cladding slab, upper bound raised to 700 K; 100 steps.
- `inputs_chk` / `inputs_restart`: the closure deck through a checkpoint
  at step 100.
- `inputs_solar`: the sun rising over the cube at Boulder on the solstice,
  from 05:00 solar time for 90 minutes; dumped every minute.

## What is checked

1. The residual of the balance on every face at every step is below 1e-3
   W/m2 (it is 1e-8; Newton takes at most 2-3 iterations).
2. The stored H, LW_net and G equal their formulas evaluated at the stored
   skin temperature (1e-8 W/m2).
3. The slab energy changes by exactly dt (G - G_bottom) every step, to the
   precision of the dump.
4. Over the run, radiation in (SW_abs + eps Q_ext + LW_net) equals the heat
   convected to the air, stored in the slab and conducted to the interior,
   to within the summed residual.
5. An independent Python model with its own Newton solve and a dense
   implicit slab, driven by the forcing dumped at every step, reproduces
   the skin temperature of every face to 1e-9 K.
6. The sunlit roof and south wall end above 302 K and the south wall is
   warmer than the north.
7. The external-flux run passes the same closure checks with every face
   past the default 380 K bound.
8. After a restart the skin and slab agree with the straight run to 1e-3 K
   (1e-5 K in practice; the immersed forcing's atmosphere is not bit-exact
   through a restart, see Phase5_Ground).
9. In the sunrise run the east wall warms before the roof, ends the warmest
   face, and the west wall stays the coolest of the walls; every face
   stays within the bounds. This check caught a mirrored solar azimuth in
   the phase 2 provider, which the noon check could not see.

## Reference output (4 ranks)

```
== closure under a fixed sun (4 ranks, 200 steps)
[IBSEB] lev=0 step=199 t=100 faces=116 (x=40 y=40 z=36) buildings=1 area=11600 m2 T_skin_min=299.9582182 T_skin_max=305.548914 SW_abs_mean=238.4418689 SW_abs_max=477.5777772 shadow_frac=0.05172413793 LW_net_mean=-59.50764834 H_mean=19.23134446 G_mean=159.7028761 Q_ext_mean=0 resid_max=9.366431186e-09 H_total_W=223083.5957
  200 dumps, steps 0-199, 116 faces
  balance residual on every face at every step: PASS (max 2.92e-08 W/m2, newton iterations at most 2)
  H, LW_net and G consistent with T_skin: PASS (max diff H 1.3e-08, LW 3.3e-09, G 8.0e-08 W/m2)
  slab energy change = dt (G - G_bottom) every step: PASS (max diff 2.0e-04 J/m2 (dump precision))
  closure: radiation in = convected + stored + conducted out: PASS (in 2.0749e+08, to air 2.2579e+07, stored -4.0224e+08, to interior 5.8714e+08 J; gap 2.84e-01, summed residual 2.57e-03)
  independent Newton + dense slab reproduces T_skin: PASS (max |dT| 1.17e-09 K over 199 steps)
  mean T_skin by orientation, first -> last dump: roof 303.39->304.15 K, east 300.34->300.39 K, west 300.34->300.39 K, north 300.28->300.34 K, south 304.40->305.21 K
  sunlit roof and south wall warm above the 300 K start, south warmer than north: PASS (roof 304.15 K, south 305.21 K, north 300.34 K)
closure: PASS
== external flux with the bound raised (4 ranks, 100 steps)
[IBSEB] lev=0 step=99 t=50 faces=116 (x=40 y=40 z=36) buildings=1 area=11600 m2 T_skin_min=381.9245916 T_skin_max=443.7037765 SW_abs_mean=0 SW_abs_max=0 shadow_frac=0 LW_net_mean=-886.9144032 H_mean=667.7245681 G_mean=1145.361029 Q_ext_mean=3000 resid_max=4.996991265e-09 H_total_W=7745604.98966
  100 dumps, steps 0-99, 116 faces
  balance residual on every face at every step: PASS (max 4.93e-08 W/m2, newton iterations at most 3)
  H, LW_net and G consistent with T_skin: PASS (max diff H 2.1e-08, LW 1.3e-08, G 1.6e-07 W/m2)
  slab energy change = dt (G - G_bottom) every step: PASS (max diff 3.0e-07 J/m2 (dump precision))
  closure: radiation in = convected + stored + conducted out: PASS (in 1.1172e+09, to air 3.5891e+08, stored 2.7431e+08, to interior 4.8400e+08 J; gap -1.04e-03, summed residual 1.49e-03)
  independent Newton + dense slab reproduces T_skin: PASS (max |dT| 9.94e-10 K over 99 steps)
  mean T_skin by orientation, first -> last dump: roof 326.13->406.03 K, east 329.34->419.64 K, west 329.34->417.32 K, north 326.14->392.53 K, south 326.14->392.09 K
  external flux on every face equals the input: PASS (3000.0 W/m2)
  every face above the default bound of 380.0 K: PASS (T_skin 381.9-443.7 K)
qext: PASS
== through a checkpoint (4 ranks)
  same faces: PASS (116 faces)
  skin and slab after the restart: PASS (max |dT_skin| 1.4e-05 K, max |dT_slab| 9.9e-07 K)
  H after the restart: PASS (max diff 1.2e-03 W/m2)
  G after the restart: PASS (max diff 1.1e-03 W/m2)
  LW_net after the restart: PASS (max diff 8.3e-05 W/m2)
  SW_abs after the restart: PASS (max diff 0.0e+00 W/m2)
restart: PASS
== sunrise over the cube (4 ranks, 10800 steps)
  step, mean T_skin east / roof / south / west / north [K]:
         0  299.76  299.63  299.72  299.70  299.74
      1200  299.99  299.55  299.65  299.62  299.80
      2400  300.69  299.70  299.72  299.67  300.12
      3600  301.76  299.98  299.84  299.78  300.58
      4800  303.10  300.35  299.97  299.91  301.11
      6000  304.61  300.81  300.13  300.06  301.65
      7200  306.23  301.33  300.28  300.22  302.17
      8400  307.89  301.96  300.43  300.37  302.67
      9600  309.57  302.60  300.57  300.51  303.05
     10680  311.06  303.22  300.69  300.64  303.36
  east wall warms before the roof: PASS (onset east 1800, roof 4200, west 7200)
  east wall ends warmest, west wall coolest of the walls: PASS (final means above)
  every face within the bounds: PASS (299.4-311.5 K)
solar: PASS
ALL PASS
```
