# SEB/Phase3_Longwave

Phase 3 of the immersed-boundary surface energy balance: the sky, ground
and building view fractions of every face by hemisphere sampling, and the
longwave balance through them. The skin temperature still does not evolve.

```
./run_longwave.sh /path/to/erf_exec        # NP=4 by default
python3 plot_longwave.py plt00002          # writes plots/*.png
```

## The scenario

The two boxes of phase 2: a tall one (80 m core at 140 m with a 5 m rim at
70 m, as the embedded-boundary reader makes it) and a short one (40 m core
at 40 m) 40 m east of it, on a 5 m grid, two steps. The sun is fixed 30 deg
above the western horizon so the shortwave of phase 2 is also exercised.
128 rays per face (16 azimuths by 8 elevations).

## The decks

- `inputs_fixed`: the sky longwave given as 300 W/m2, ground at 300 K with
  emissivity 0.95, faces with emissivity 0.9.
- `inputs_gray`: the sky longwave as 0.83 sigma T_air^4 with the air
  temperature of each face's fluid cell.

## What is checked

1. The three view fractions of every face sum to one.
2. They equal an independent Python hemisphere sampling on every face: the
   same stratified directions and an independently written ray walk over
   column tops rebuilt from the roof faces of the dump.
3. Analytic values: the tall core roof sees only sky; no roof sees the
   ground; the tall west wall above its rim, facing open ground, sees
   exactly half sky, and the 5 m ledge below it fills a quarter of the view
   of the rows just above it but little from 50 m higher (the few percent
   that remain there are rays wrapping round the periodic domain and
   meeting the box's own east wall from behind); the short box's core west
   wall sees the tall box in most of its view.
4. The longwave formulas on every face, with the sky term fixed or gray,
   and that the air temperature read from the fluid cells is plausible.
5. The fixed deck on one rank equals the four-rank dump.

## Reference

```
== fixed (4 ranks)
[IBSEB DEBUG] lev=0 view fractions: 128 rays per face, mean f_sky=0.530909547 f_ground=0.3711893157 f_bldg=0.09790113723
  f_sky + f_ground + f_bldg = 1 on every face              PASS 
  fractions = independent hemisphere sampling, every face  PASS mismatches 0/2616
  tall core roof: f_sky = 1                                PASS n=256
  every roof: f_ground = 0                                 PASS 
  tall west wall above the rim: f_sky = 1/2                PASS f_sky 0.500-0.500
  tall west wall: ledge in view fades with height above it PASS f_bldg 0.062 at z > 120 m, 0.227 just above the ledge
  short core west wall sees the tall box: f_bldg > 0.2     PASS f_bldg 0.625-0.812
  LW_in = f_sky LW_sky + f_ground eps_g sigma Tg^4 + f_bldg sigma Ts^4 PASS 
  LW_net = eps (LW_in - sigma Ts^4)                        PASS 
== gray (4 ranks)
[IBSEB DEBUG] lev=0 view fractions: 128 rays per face, mean f_sky=0.530909547 f_ground=0.3711893157 f_bldg=0.09790113723
  f_sky + f_ground + f_bldg = 1 on every face              PASS 
  fractions = independent hemisphere sampling, every face  PASS mismatches 0/2616
  tall core roof: f_sky = 1                                PASS n=256
  every roof: f_ground = 0                                 PASS 
  tall west wall above the rim: f_sky = 1/2                PASS f_sky 0.500-0.500
  tall west wall: ledge in view fades with height above it PASS f_bldg 0.062 at z > 120 m, 0.227 just above the ledge
  short core west wall sees the tall box: f_bldg > 0.2     PASS f_bldg 0.625-0.812
  LW_in = f_sky LW_sky + f_ground eps_g sigma Tg^4 + f_bldg sigma Ts^4 PASS 
  LW_net = eps (LW_in - sigma Ts^4)                        PASS 
  gray sky: T_air read from the fluid cells is plausible (250-320 K) PASS 289.3-290.6 K
== fixed on 1 rank against 4 ranks
  rank independence of the face dump: PASS (2616 vs 2616 faces)
ALL PASS
```
