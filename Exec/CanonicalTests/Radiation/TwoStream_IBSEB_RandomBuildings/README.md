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
- **The faces' sky.** `erf.ibseb.sun_mode = two_stream`: the faces take the sun position from the
  two-stream radiation (its orbital formula has no equation of time, so the prescribed `solar` sun
  would sit 1.2° off it on this date). `erf.ibseb.radiation = two_stream`: they take their
  radiation from the two-stream columns, each face from the column of its fluid cell at its own
  height: the direct beam (as a direct normal irradiance, the beam over cos z), the diffuse sky
  and the longwave down there, and the ground's shortwave and longwave up, which carry the
  force-restore ground's albedo, emissivity and skin temperature. The column is clear and purely
  absorbing (0.002 per 10 m layer, no scattering), so the beam at height z is
  S0 cos z exp(-0.002 (100 - z / 10 m) / cos z) and there is no diffuse sky.

## Running it

```
python3 make_buildings.py                # the files are in the directory already; rerun for another set
mpirun -np 4 erf_exec inputs             # about 2 h 30 min on 4 ranks of a laptop
python3 check_random_buildings.py        # exit code = number of failed checks
python3 plot_random_buildings.py         # random_buildings_{map,tskin,theta}.png
```

## What is checked

`check_random_buildings.py` reads `buildings.csv`, the per-building report `ibseb_buildings.csv`,
the face dumps (`faces/set.*` on level 0, `faces/set.lev1.*` on level 1), `radiation_diag.csv` and
the last 2D plotfile:

1. Both levels hold all 8 buildings, numbered alike, their roofs centred on `buildings.csv`.
2. The balance closes on every building and level (residual below 1e-3 W/m²).
3. The faces see the two-stream sun and the columns' beam: on every sunlit, unshadowed roof of
   level 0 at every face dump, the top-of-atmosphere irradiance its direct shortwave implies at
   its own height equals the sweep's, to 1e-5, with no diffuse sky and the sun east of the
   meridian. A report row and a dump carry the time at the end of their step and the sun the step
   used, from its start, so they are compared with the sweep one step earlier.
4. The sun climbs and the shadows shorten on both levels, from the first report after the initial
   one (which comes before any sweep and carries no sun on the faces).
5. Every building warms on both levels.
6. The two levels agree on what they resolve alike: each building's top roof within 0.5 K and its
   walls within 1 K at the last face dump. The whole-building means are reported, not checked (see
   below).
7. The open ground (away from the footprints) warms on both levels.

It also prints, for every building and level, the area, absorbed shortwave and skin temperature
of its top roof, its ledges and its walls; the reading below rests on those numbers.

## Results

On 4 MPI ranks of a laptop shared with other runs, 10,800 level-0 steps in about 2 h 30 min of
run time. The balance costs 0.25 ms per step on level 0 and 0.11 ms per step on level 1 on the
slowest rank.

All 7 checks pass:

| Check | Result |
|---|---|
| 1. Buildings | 8 and 8, roof centres on `buildings.csv` |
| 2. Residual | largest 3.3e-8 W/m² over 2,896 rows |
| 3. The sun and the beam | 11,034 sunlit roof faces at 5 heights over the dumps, largest difference 4.8e-7; no diffuse sky; azimuth 99.5° to 156.6° |
| 4. The sun climbs | zenith 52.7° to 24.8°; shadowed face area 5.3 % to 3.2 % on level 0, 6.5 % to 4.2 % on level 1 |
| 5. Warming | final means 314.0–321.7 K on level 0, 312.6–319.6 K on level 1 (from 300 K) |
| 6. Levels | level 1 minus level 0: top roofs -0.12 to -0.01 K, walls -0.64 to +0.33 K (whole buildings -2.16 to -0.24 K) |
| 7. Open ground | 310.52 K in the refined area (level 1), 310.42 K outside it (level 0) at 11:20 |

The sun of the faces, with the mean over the level-0 faces of the direct normal irradiance each
takes at its height (the first report with a sun is one minute in):

| UTC | Zenith | Azimuth | Direct normal | Diffuse, horizontal |
|---|---|---|---|---|
| 15:01 | 52.5° | 99.7° | 958 W/m² | 0 |
| 16:00 | 41.6° | 112.1° | 1017 W/m² | 0 |
| 17:00 | 31.8° | 129.8° | 1050 W/m² | 0 |
| 18:00 | 24.8° | 156.6° | 1065 W/m² | 0 |

The buildings at 18:00 UTC, level 0 / level 1:

| Building | Footprint | Height | Material | Faces | Face area [m²] | Mean skin [K] | Absorbed SW [W/m²] | Net LW [W/m²] | H [W/m²] |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 20 m × 40 m | 60 m | timber | 52 / 120 | 12000 / 12000 | 318.1 / 316.4 | 225 / 182 | -76 / -64 | 110 / 89 |
| 2 | 40 m × 20 m | 30 m | timber | 34 / 68 | 8400 / 6800 | 321.7 / 319.6 | 306 / 237 | -94 / -78 | 159 / 110 |
| 3 | 20 m × 20 m | 60 m | brick | 37 / 88 | 8400 / 8800 | 314.0 / 312.6 | 243 / 207 | -50 / -40 | 72 / 57 |
| 4 | 60 m × 20 m | 20 m | brick | 35 / 68 | 9200 / 6800 | 317.3 / 316.8 | 408 / 378 | -84 / -72 | 174 / 153 |
| 5 | 40 m × 60 m | 10 m | concrete | 34 / 72 | 10000 / 7200 | 317.7 / 317.5 | 525 / 509 | -95 / -90 | 187 / 183 |
| 6 | 60 m × 20 m | 50 m | timber | 59 / 132 | 14000 / 13200 | 318.3 / 318.1 | 230 / 221 | -76 / -73 | 117 / 111 |
| 7 | 20 m × 20 m | 50 m | timber | 33 / 76 | 7600 / 7600 | 318.5 / 317.9 | 219 / 204 | -82 / -74 | 96 / 87 |
| 8 | 20 m × 20 m | 20 m | brick | 21 / 36 | 5200 / 3600 | 316.8 / 316.3 | 356 / 330 | -77 / -65 | 139 / 125 |

The parts of each building at 18:00 UTC, level 0 / level 1 (from the checker's breakdown):

| Building | Top roof: SW [W/m²], skin [K] | Ledges: area [m²], SW [W/m²], skin [K] | Walls: SW [W/m²], skin [K] |
|---|---|---|---|
| 1 | 586 / 586, 324.0 / 324.0 | 2400 / 1600, 483 / 327, 321.9 / 315.9 | 122 / 124, 316.5 / 315.9 |
| 2 | 582 / 582, 324.7 / 324.7 | 2400 / 1200, 579 / 387, 326.8 / 318.3 | 137 / 142, 318.9 / 319.0 |
| 3 | 683 / 683, 323.1 / 323.0 | 1600 / 1200, 507 / 395, 319.6 / 315.4 | 149 / 149, 312.0 / 311.6 |
| 4 | 677 / 677, 323.3 / 323.3 | 3200 / 1600, 676 / 676, 323.9 / 324.3 | 163 / 170, 311.5 / 311.8 |
| 5 | 724 / 724, 321.8 / 321.7 | none | 172 / 172, 310.5 / 310.8 |
| 6 | 584 / 584, 324.1 / 324.0 | 3200 / 2000, 362 / 378, 317.5 / 317.7 | 142 / 146, 317.9 / 317.5 |
| 7 | 584 / 584, 325.5 / 325.5 | 1600 / 1200, 434 / 435, 321.3 / 320.7 | 132 / 133, 317.2 / 316.8 |
| 8 | 677 / 677, 323.4 / 323.3 | 1600 / 800, 676 / 676, 324.2 / 324.1 | 155 / 156, 312.3 / 312.6 |

What the numbers show:

- **Roofs absorb four to five times what walls do.** The top roofs absorb 582–724 W/m² at 18:00
  UTC and the walls 122–172 W/m². So the buildings whose faces are mostly roof absorb the most per
  square metre: the 10 m concrete block (5), whose faces are 64 % roof on level 0, absorbs
  525 W/m² on average; the 60 m towers (1, 3), mostly wall, 182–243 W/m².
- **A roof absorbs by material and height, the same on both levels.** A top roof's absorbed
  shortwave is set by its albedo, the sun and the air above it: timber 582–586 W/m², the higher
  roofs a little more (less absorbing air above); brick 677–683 W/m²; concrete 724 W/m².
- **The two levels resolve the top roofs and the walls alike.** Their skin temperatures agree to
  0.12 K on the top roofs and 0.64 K on the walls.
- **The levels differ at the buildings' edges.** Each level steps a building's edge down over one
  of its cells, as upward-facing ledges. On level 0 that is one ledge 10 m up and 20 m wide; on
  level 1, ledges 10 m wide, 10 m to 30 m up. On buildings 1, 2 and 3, level 1's ledges absorb
  327–395 W/m² against 483–579 W/m² on level 0, and run 4–9 K cooler. On the others both levels'
  ledges absorb within 16 W/m² of each other, and on the 20 m buildings (4, 8) they are fully
  sunlit. The 10 m block (5) has no ledge: its step is at roof height and counts as roof. The
  ledges cover more of the face area on level 0 (19–35 % against 13–24 % on level 1), and absorb
  more than the walls. With the lower absorption on buildings 1–3, this sets the whole-building
  differences: a lower mean absorbed shortwave on level 1 for every building, and a skin
  0.2–2.2 K cooler.
- **The longwave from the columns warms the buildings, the walls most.** Against the run of this
  case before the two-stream provider (the same beam through matched clear-sky settings, a gray
  sky of emissivity 0.83 and a fixed 300 K ground), the absorbed shortwave agrees within 1 % on
  every building and level 0, while the buildings end 1.2–4.1 K warmer and their walls 1.8–5.2 K
  warmer. The faces now take the column's sky and the force-restore ground (about 310 K by 11:20
  local) through their ground view; that run does not separate the two.
- **The open ground is the same on both levels.** It reaches 310.52 K in the refined area and
  310.42 K outside it. Under and next to the footprints the ground balance runs hotter (312.09 K on
  level 1): the two-stream columns do not stop at the buildings, so it still takes the sun, and it
  loses less sensible heat there (236 W/m² against 330 W/m² over the open ground).

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
- The two-stream diagnostics CSV leaves `T_s_mean` empty (`nan`) with the prognostic balance; the
  ground temperatures above are from the 2D plotfiles.
