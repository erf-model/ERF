# SEB/Phase2_Shortwave

Phase 2 of the immersed-boundary surface energy balance: shortwave on the
building faces with ray-cast shadowing, from the prescribed radiation
provider. The skin temperature still does not evolve.

```
./run_shortwave.sh /path/to/erf_exec        # NP=4 by default
python3 plot_shortwave.py plt00002          # writes plots/*.png
```

## The scenario

A 640 m cube domain on a 5 m grid with two immersed-forcing boxes from
`two_buildings_5m_128x128.txt` (made by `make_buildings.py`): a tall one,
80 x 80 x 140 m at x = 280-360 m, and a short one, 40 x 40 x 40 m at
x = 400-440 m, both centred on y = 320 m, so the gap between the tall box's
east wall and the short box's west wall is 40 m. Two steps; the face dump at
step 2 is what the checker reads.

## The decks

- `inputs_zen60`: fixed sun in the west, 30 deg above the horizon, direct
  normal 800 W/m2, diffuse 100 W/m2. The tall box's shadow reaches
  140 - 40 tan(30) = 117 m up the short box's west wall, which is 40 m tall,
  so that wall is entirely in shadow.
- `inputs_zen20`: the sun 70 deg above the horizon; the shadow reaches 30 m,
  so the six lowest rows of the short west wall are shadowed and the top two
  and the roof are lit.
- `inputs_solar`: sun from the site and time, Boulder at solar noon on the
  June solstice, where the zenith must be the latitude minus the declination
  and the sun due south.

## What the geometry really is

The embedded-boundary reader interpolates the nodal heightmap, so each
building edge becomes a one-cell step: the tall box is an 80 m full-height
core at 140 m with a 5 m rim at 70 m and 5 m corner stubs, and the short box
a 40 m core at 40 m with a rim at 20 m. The checker does not assume plane
walls; it reads the geometry from the face dump.

## What is checked

1. **Shadow flag against an independent ray cast.** The column tops are
   rebuilt from the roof faces of the dump and a ray is cast from every
   face toward the sun in Python; the code's flag must match on every face
   (2616 faces, 644 in shadow at zen20).
2. **Fluxes on every face.** Direct = DNI max(0, n.s) (1 - shadow); the
   placeholder diffuse; absorbed = (1 - albedo) times the sum.
3. **Analytic spot checks.** The tall core roof and the tall west wall are
   unshadowed with the exact incidence; on the short box's core west wall
   the shadow stops at H - gap tan(elevation) with the 40 m gap from the
   tall core's east edge (30.1 m at zen20: the rows at 22.5 and 27.5 m are
   dark, 32.5 and 37.5 m lit). The rim wall 5 m nearer is fully shadowed,
   and on the tall roof only the east rim, stepped below the core, is.
4. The zen60 dump on one rank equals the four-rank dump; the solar deck's
   zenith and azimuth match the solstice-noon formulas.

## Reference

```
== zen60 (4 ranks)
[IBSEB DEBUG] lev=0 sun: zenith=60 deg azimuth=270 deg s=(-0.8660254038,-1.590862858e-16,0.5) DNI=800 W/m2 diffuse_h=100 W/m2
  shadow flag = independent ray cast, every face       PASS mismatches 0/2616, shadowed 738
  direct = DNI max(0, n.s) (1 - shadow)                PASS 
  diffuse = f_sky D + f_ground a_g (DNI cos z + D)     PASS 
  absorbed = (1 - albedo) (direct + diffuse)           PASS 
  tall core roof: unshadowed, direct = DNI cos z       PASS n=256
  tall west wall: unshadowed, direct = DNI sin z       PASS n=476
  short core west wall: shadow where z < 116.9 m (gap 40 m) PASS shadowed 38/38
== zen20 (4 ranks)
[IBSEB DEBUG] lev=0 sun: zenith=20 deg azimuth=270 deg s=(-0.3420201433,-6.282808107e-17,0.9396926208) DNI=800 W/m2 diffuse_h=100 W/m2
  shadow flag = independent ray cast, every face       PASS mismatches 0/2616, shadowed 644
  direct = DNI max(0, n.s) (1 - shadow)                PASS 
  diffuse = f_sky D + f_ground a_g (DNI cos z + D)     PASS 
  absorbed = (1 - albedo) (direct + diffuse)           PASS 
  tall core roof: unshadowed, direct = DNI cos z       PASS n=256
  tall west wall: unshadowed, direct = DNI sin z       PASS n=476
  short core west wall: shadow where z < 30.1 m (gap 40 m) PASS shadowed 22/38
== solar (4 ranks)
[IBSEB DEBUG] lev=0 sun: zenith=16.59919815 deg azimuth=181.0445156 deg s=(-0.00520763337,-0.2856274863,0.9583265725) DNI=907.4975233 W/m2 diffuse_h=1
  solar noon Boulder: zenith 16.60 (expected 16.60), azimuth 181.0 (expected 180) -> PASS
  clear-sky beam on roofs: max direct 870 W/m2 (> 500) -> PASS
== zen60 on 1 rank against 4 ranks
  rank independence of the face dump: PASS (2616 vs 2616 faces)
ALL PASS
```
