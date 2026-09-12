# SEB/Phase7_IsolatedBuilding

Canonical case for the immersed-boundary surface energy balance: one
building through a full day. It exercises everything the eight phases
built (shadow and incidence, view fractions and gray sky, the wall
function, the slab, the prognostic skin) on the simplest geometry, with
checks on the sequence of the day and on two integrals, so a change to
any part of the balance shows up here.

```
./run_isolated.sh /path/to/erf_exec                      # NP=4 by default, about two hours
python3 plot_isolated.py ibseb_day.csv faces/day --plotfile plt100800    # writes plots/*.png
```

## The scenario

A 40 m cube, an exact box on cell faces (`eb2.geometry = box`), of 30 cm concrete (k = 1.5 W/m/K, rho c_p = 2e6 J/m3/K, ten
layers, albedo 0.25, emissivity 0.92, 295 K inside, faces and slabs
starting at the air temperature) at the centre of a
320 m periodic domain on a 10 m grid, in a 3 m/s westerly, neutral at
300 K, with a MOST ground (z0 = 0.1 m) held at 300 K. Boulder on the June
solstice, from midnight solar time for 24 hours: clear sky from the
prescribed provider (Bird direct beam with transmission 0.7, Liu-Jordan
diffuse), gray sky longwave with emissivity 0.83, and a 300 K ground seen
by the walls. The domain is small so the day runs in about two hours on four ranks
(130 minutes on a laptop); the balance itself costs nothing at this size.

Output: `ibseb_day.csv` with the building means, the skin temperature
range, the largest residual and the sun every minute; `faces/day.step*`
with every face every minute; plotfiles every two hours; a checkpoint at
the end.

## What happens through the day

Mean skin temperature by orientation, and the air next to the roof, in K:

```
hour   air    roof   east   south  west   north
 0.0  299.6  299.3  299.6  299.6  299.6  299.6
 4.0  301.5  296.5  297.2  297.7  297.4  297.7
 6.0  301.9  299.6  303.4  298.7  298.5  300.5
 8.0  302.2  311.6  318.7  301.0  300.8  302.8
10.0  302.6  325.7  322.8  307.2  303.0  303.9
12.0  303.0  336.1  315.9  313.1  304.6  305.1
14.0  303.1  339.7  312.1  314.4  314.6  305.6
16.0  303.3  335.0  309.7  310.0  324.2  305.1
18.0  303.3  323.9  306.8  306.3  322.2  306.2
20.0  303.3  312.4  302.7  302.2  309.2  301.1
22.0  303.3  307.1  300.5  300.1  304.7  299.1
```

- **Night.** Every face radiates to a gray sky at 0.83 emissivity and ends
  the night below the air; the roof, which sees the whole sky, is the
  coldest at 04:30, 5.3 K below the air and 1 K below the walls. The
  interior at 295 K draws a little heat as well (G of -25 to -35 W/m2).
- **Sunrise (04:35 solar time).** The east wall is the first face to
  rise above the air, at 05:51, forty minutes before the roof; the
  north wall follows because the solstice sun rises north of east, then
  the south wall at 08:40 and the west wall at 09:43.
- **Morning.** The east wall peaks at 09:39 at 323 K and then falls as
  the beam leaves it. The roof climbs all morning, storing 400 to 480 W/m2
  in the concrete (G) while the net longwave loss grows to -300 W/m2.
- **Afternoon.** The roof peaks at 13:55 at 340 K, two hours after the
  sun, with the absorbed beam already falling from its 800 W/m2 noon
  value; the south wall peaks at 13:23 at 315 K and is 9 K warmer than
  the north wall at 13:00; the west wall peaks at 16:45 at 325 K.
- **Evening.** Conduction turns around at 18:10: the slab, which stored
  16.5 MJ/m2 through the day on the roof, gives it back and keeps the skin
  above the air until after 23:00 on the roof and the west wall.
- **The atmosphere** stays clean on the box: fluid theta within 0.4 K of
  neutral at 12:00, vertical velocity under 0.4 m/s, a warm plume of 0.3 K
  rising downwind of the cube at 14:00 (`theta_slice.png`).
- **A caveat.** The sensible flux off the 340 K roof is only 10 to 30
  W/m2 because the wall function is neutral with a fixed roughness and the
  roof wind is light in the separation zone; the roof is hot mostly
  because it cannot convect. The stability functions and a convective
  velocity scale at low wind are the first items of the next PR (see the
  plan file); the day will be rerun with them on.

## What is checked

1. The balance residual stays below 1e-3 W/m2 all day.
2. At night every face loses longwave, and before sunrise the roof,
   which sees the whole sky, is colder than the air and than the walls.
3. After sunrise the east wall is the first face to rise above the air
   next to it (every face is colder than the air at dawn).
4. The peaks come in the order east wall, roof, west wall, with the roof
   peaking in the early afternoon (12:00 to 15:30 solar time).
5. At 13:00 the south wall is warmer than the north wall.
6. The shortwave absorbed by the core roof over the day equals the
   clear-sky formulas integrated independently in Python at the dump
   times, to 1 percent.
7. The slab energy of every face changes over the day by the integrated
   conduction, to 3 percent of the integrated |G| (trapezoid on the
   minute dumps).

## Not in this case, and why the box

The ground has no balance of its own: the development branch's
land-surface models have no radiation, so the ground stays at 300 K for
the walls' longwave term and supplies drag only. A 1 s step trips the
dycore's floating-point check on this grid with the balance off as well,
so the day runs at 0.5 s.

The cube is an exact box rather than the nodal height map of the earlier
phases: the height map's one-cell ramp leaves corner cells with 1 to 20
percent solid, and on those the immersed forcing develops a vertical
checkerboard in temperature after about two hours of this run, with the
balance off as well (see the findings in
`Source/ImmersedBoundarySEB/IBSEB_DEVELOPMENT.md`). The box has no partial
cells and 80 plane faces.

## Reference output (4 ranks)

```
== isolated building, 24 h (4 ranks)
  balance residual all day: PASS (max 3.7e-08 W/m2 over 1441 rows)
  hour   T_air  roof   east   south  west   north  [K]
    0.0  299.56 299.32 299.56 299.58 299.56 299.58
    2.0  300.81 297.35 298.04 298.31 298.14 298.40
    4.0  301.48 296.51 297.24 297.68 297.42 297.70
    6.0  301.86 299.55 303.41 298.74 298.48 300.50
    8.0  302.21 311.57 318.68 301.01 300.84 302.82
   10.0  302.60 325.72 322.75 307.16 302.98 303.88
   12.0  302.95 336.12 315.94 313.08 304.55 305.11
   14.0  303.14 339.67 312.07 314.35 314.59 305.60
   16.0  303.29 335.03 309.69 309.95 324.19 305.08
   18.0  303.26 323.87 306.83 306.34 322.17 306.21
   20.0  303.28 312.37 302.71 302.22 309.19 301.14
   22.0  303.32 307.05 300.46 300.09 304.67 299.11
  night: every face loses longwave; before sunrise the roof is colder than the air and than the walls: PASS (max LW_net at night -25.4 W/m2, at 04:30 roof - air -5.33 K, roof - walls -1.07 K)
  east wall is the first face to rise above the air after sunrise: PASS (east 5.85 h, roof 6.50 h, north 6.70 h, south 8.67 h, west 9.72 h)
  peaks in the order east, roof, west with the roof in the early afternoon: PASS (east 9.65 h (322.9 K), roof 13.92 h (339.7 K), west 16.75 h (325.5 K), south 13.38 h (314.6 K), north 17.68 h (306.3 K))
  south wall warmer than the north wall at 13:00: PASS (south 314.47 K, north 305.46 K)
  daily shortwave absorbed by the core roof vs the clear-sky formulas integrated in Python: PASS (ERF 24.882 MJ/m2, Python 24.882 MJ/m2)
  slab energy over the day vs the integrated conduction (trapezoid on minute dumps): PASS (worst face 0.01 % of the integrated |G| (16.54 MJ/m2))
isolated building: PASS
ALL PASS
```
