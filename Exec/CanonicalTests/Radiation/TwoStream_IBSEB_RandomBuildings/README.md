# Two-stream radiation over random buildings, on two levels

The third of the two-stream test cases, after the Noah-MP and force-restore comparison in
`../TwoStream_NoahMP_vs_ForceRestore`. A seeded random set of buildings stands on the ground of
that case: the two-stream radiation heats the atmosphere and drives the force-restore ground, and
the immersed-boundary surface energy balance (`erf.ibseb`) heats the walls and roofs. The faces see
the same sun as the two-stream columns, and the buildings sit on a refined level.

## Setup

- **Buildings.** `make_buildings.py` places 8 box buildings by seeded rejection sampling (seed 7)
  in x, y = 200–440 m: footprints of 20–60 m on 20 m corners, heights of 10–60 m on 10 m steps,
  streets at least 60 m wide, and a material (concrete, brick or timber from `materials.csv`) each.
  It writes the height map `buildings.txt`, the table `buildings.csv` and `inputs_buildings`, which
  the deck reads through `FILE =`. Run it with another `--seed` or `--n` for another set.
- **Grid.** 640 m × 640 m × 1 km, periodic in x and y. Level 0 is 20 m × 20 m × 10 m; level 1 is
  10 m × 10 m × 10 m over x, y = 160–480 m and spans the depth. It holds every building with at
  least one level-0 cell around it; `python3 ../../SEB/ibseb_refinement_box.py inputs` confirms it.
- **Step.** 1 s on level 0, two 0.5 s steps on level 1.
- **Atmosphere.** Dry, 300 K up to 300 m and 10 K/km above, a 5 m/s geostrophic westerly at 40° N,
  Smagorinsky LES, 0.1 K random perturbations.
- **Sun.** 2024-08-05 from 15:00 to 18:00 UTC at 40° N, 100° W, about 08:20 to 11:20 local solar
  time.
- **Ground.** The two-stream radiation with the force-restore balance of case 2: two-way with the
  surface layer, albedo 0.21, emissivity 0.994, deep temperature 290 K. The balance covers every
  column, under the buildings too, so the checks leave the footprints out.
- **Buildings' balance.** Prognostic skin, 10 slab layers, the convective velocity scale and the
  stability functions on, 16 × 8 rays for the view fractions.
- **The faces' sky.** `erf.ibseb.sun_mode = two_stream`: the faces take the sun position and the
  top-of-atmosphere irradiance S0 from the two-stream radiation. The radiation's orbital formula
  has no equation of time, so the prescribed `solar` sun would sit 1.2° off it on this date. The
  two-stream column is clear and purely absorbing with an optical depth of 0.2, so the faces' beam
  crosses the same atmosphere: `sw_transmission = exp(-0.2)` and `sw_diffuse_coeff = 0`, which
  give a direct normal irradiance of S0 exp(-0.2 / cos z) and no diffuse sky.

## Running it

```
python3 make_buildings.py                # the files are in the directory already; rerun for another set
mpirun -np 4 erf_exec inputs             # about 2 h 15 min on 4 ranks of a laptop
python3 check_random_buildings.py        # exit code = number of failed checks
python3 plot_random_buildings.py         # random_buildings_{map,tskin,theta}.png
```

## What is checked

`check_random_buildings.py` reads `buildings.csv`, the per-building report `ibseb_buildings.csv`,
the face dumps (`faces/set.*` on level 0, `faces/set.lev1.*` on level 1), `radiation_diag.csv` and
the last 2D plotfile:

1. Both levels hold all 8 buildings, numbered alike, their roofs centred on `buildings.csv`.
2. The balance closes on every building and level (residual below 1e-3 W/m²).
3. The faces see the two-stream sun: the top-of-atmosphere irradiance their direct beam implies
   equals the sweep's, to 1e-5, with the sun east of the meridian. A report row carries the time at
   the end of its step and the sun the step used, from its start, so it is compared with the sweep
   one step earlier.
4. The sun climbs and the shadows shorten on both levels.
5. Every building warms on both levels.
6. The two levels agree on what they resolve alike: each building's top roof within 0.5 K and its
   walls within 1 K at the last face dump. The whole-building means are reported, not checked (see
   below).
7. The open ground (away from the footprints) warms on both levels.

It also prints, for every building and level, the area, absorbed shortwave and skin temperature
of its top roof, its ledges and its walls; the reading below rests on those numbers.

## Results

On 4 MPI ranks of a laptop, 10,800 level-0 steps. The balance costs 0.36 ms per step on level 0
and 0.12 ms per step on level 1 on the slowest rank.

All 7 checks pass:

| Check | Result |
|---|---|
| 1. Buildings | 8 and 8, roof centres on `buildings.csv` |
| 2. Residual | largest 3.3e-8 W/m² over 2,896 rows |
| 3. The sun | largest difference 4.8e-7 over 180 times; azimuth 99.5° to 156.6° |
| 4. The sun climbs | zenith 52.7° to 24.8°; shadowed face area 5.3 % to 3.2 % on level 0, 6.5 % to 4.2 % on level 1 |
| 5. Warming | final means 311.9–317.6 K on level 0, 310.5–316.3 K on level 1 (from 300 K) |
| 6. Levels | level 1 minus level 0: top roofs -0.08 to +0.11 K, walls -0.63 to +0.25 K (whole buildings -2.97 to -0.26 K) |
| 7. Open ground | 310.53 K in the refined area (level 1), 310.47 K outside it (level 0) at 11:20 |

The sun of the faces at the start of each hour:

| UTC | Zenith | Azimuth | Direct normal | Diffuse, horizontal |
|---|---|---|---|---|
| 15:00 | 52.7° | 99.5° | 951 W/m² | 0 |
| 16:00 | 41.6° | 112.1° | 1012 W/m² | 0 |
| 17:00 | 31.8° | 129.8° | 1045 W/m² | 0 |
| 18:00 | 24.8° | 156.6° | 1061 W/m² | 0 |

The buildings at 18:00 UTC, level 0 / level 1:

| Building | Footprint | Height | Material | Faces | Face area [m²] | Mean skin [K] | Absorbed SW [W/m²] | H [W/m²] |
|---|---|---|---|---|---|---|---|---|
| 1 | 20 m × 40 m | 60 m | timber | 52 / 120 | 12000 / 12000 | 315.0 / 313.2 | 224 / 180 | 88 / 63 |
| 2 | 40 m × 20 m | 30 m | timber | 34 / 68 | 8400 / 6800 | 317.6 / 314.6 | 305 / 236 | 173 / 127 |
| 3 | 20 m × 20 m | 60 m | brick | 37 / 88 | 8400 / 8800 | 311.9 / 310.5 | 242 / 205 | 67 / 52 |
| 4 | 60 m × 20 m | 20 m | brick | 35 / 68 | 9200 / 6800 | 315.9 / 315.1 | 407 / 377 | 150 / 134 |
| 5 | 40 m × 60 m | 10 m | concrete | 34 / 72 | 10000 / 7200 | 316.5 / 316.3 | 524 / 508 | 156 / 148 |
| 6 | 60 m × 20 m | 50 m | timber | 59 / 132 | 14000 / 13200 | 315.0 / 314.5 | 229 / 220 | 98 / 93 |
| 7 | 20 m × 20 m | 50 m | timber | 33 / 76 | 7600 / 7600 | 314.7 / 313.8 | 219 / 203 | 91 / 85 |
| 8 | 20 m × 20 m | 20 m | brick | 21 / 36 | 5200 / 3600 | 315.0 / 314.3 | 355 / 329 | 118 / 103 |

The parts of each building at 18:00 UTC, level 0 / level 1 (from the checker's breakdown):

| Building | Top roof: SW [W/m²], skin [K] | Ledges: area [m²], SW [W/m²], skin [K] | Walls: SW [W/m²], skin [K] |
|---|---|---|---|
| 1 | 578 / 578, 323.3 / 323.4 | 2400 / 1600, 482 / 325, 320.7 / 313.9 | 121 / 123, 312.6 / 312.3 |
| 2 | 578 / 578, 321.4 / 321.5 | 2400 / 1200, 578 / 385, 324.7 / 315.5 | 137 / 142, 313.7 / 313.2 |
| 3 | 674 / 674, 321.2 / 321.1 | 1600 / 1200, 506 / 393, 318.6 / 314.2 | 149 / 148, 309.6 / 309.3 |
| 4 | 674 / 674, 321.9 / 321.9 | 3200 / 1600, 674 / 674, 323.3 / 323.6 | 163 / 169, 309.5 / 309.7 |
| 5 | 722 / 722, 320.9 / 320.9 | none | 172 / 172, 308.7 / 308.9 |
| 6 | 578 / 578, 323.2 / 323.2 | 3200 / 2000, 361 / 376, 315.6 / 315.8 | 141 / 146, 313.8 / 313.2 |
| 7 | 578 / 578, 324.4 / 324.4 | 1600 / 1200, 433 / 433, 319.7 / 319.1 | 131 / 132, 312.5 / 312.0 |
| 8 | 674 / 674, 322.2 / 322.1 | 1600 / 800, 674 / 674, 323.2 / 323.2 | 155 / 156, 310.0 / 310.1 |

What the numbers show:

- **Roofs absorb four to five times what walls do.** The top roofs absorb 578–722 W/m² at 18:00
  UTC and the walls 121–172 W/m². So the buildings whose faces are mostly roof absorb the most per
  square metre: the 10 m concrete block (5), whose faces are 64 % roof on level 0, absorbs
  524 W/m² on average; the 60 m towers (1, 3), mostly wall, 180–242 W/m².
- **The roof absorbs by material, the same on both levels.** A top roof's absorbed shortwave is
  set by its albedo and the sun: 578 W/m² for timber, 674 W/m² for brick, 722 W/m² for concrete.
- **The two levels resolve the top roofs and the walls alike.** Their skin temperatures agree to
  0.11 K on the top roofs and 0.63 K on the walls.
- **The levels differ at the buildings' edges.** Each level steps a building's edge down over one
  of its cells, as upward-facing ledges. On level 0 that is one ledge 10 m up and 20 m wide; on
  level 1, ledges 10 m wide, 10 m to 30 m up. On buildings 1, 2 and 3, level 1's ledges absorb
  325–393 W/m² against 482–578 W/m² on level 0, and run 4–9 K cooler. On the others both levels'
  ledges absorb within 15 W/m² of each other, and on the 20 m buildings (4, 8) they are fully
  sunlit. The 10 m block (5) has no ledge: its step is at roof height and counts as roof. The
  ledges cover more of the face area on level 0 (19–35 % against 13–24 % on level 1), and absorb
  more than the walls. With the lower absorption on buildings 1–3, this sets the whole-building
  differences: a lower mean absorbed shortwave on level 1 for every building, and a skin
  0.3–3.0 K cooler.
- **The open ground is the same on both levels.** It reaches 310.53 K in the refined area and
  310.47 K outside it. Under and next to the footprints the ground balance runs hotter (312.12 K on
  level 1): the two-stream columns do not stop at the buildings, so it still takes the sun, and it
  loses less sensible heat there (221 W/m² against 284 W/m² over the open ground).

## Plots

`plot_random_buildings.py` draws:

- **`random_buildings_map.png`:** the last 2D plotfile from above: the open ground's skin, the
  footprints and the level-1 roofs.
- **`random_buildings_tskin.png`:** the skin temperature of every building on both levels against
  time, with the open ground.
- **`random_buildings_theta.png`:** the potential temperature 15 m up on level 1, with the wind.

## Notes

- The mixed-layer depth of the convective velocity scale (`erf.ibseb.z_i_mode = bulk_ri`) comes
  from level 0's profile on both levels.
- The faces' ground longwave takes the constant `erf.ibseb.T_ground = 300 K`, not the ground
  balance's skin, which warms to about 310 K.
- The two-stream diagnostics CSV leaves `T_s_mean` empty (`nan`) with the prognostic balance; the
  ground temperatures above are from the 2D plotfiles.
