# Two-stream radiation over random buildings, on two levels

The third of the two-stream test cases. It takes the ground of case 2
(`../TwoStream_NoahMP_vs_ForceRestore`, the two-stream radiation with the force-restore surface
balance) and stands a seeded random set of buildings on it. The two-stream radiation heats the air
and the ground. The immersed-boundary surface energy balance (`erf.ibseb`) heats the walls and roofs,
and its faces take their radiation from the same two-stream columns. The buildings sit on a
refined level.

## Setup

| Item | Setting |
|---|---|
| Buildings | 8 boxes from `make_buildings.py` (seed 7) in x, y = 200–440 m: footprints 20–60 m, heights 10–60 m, streets at least 60 m wide; concrete, brick or timber (`materials.csv`) |
| Grid | 640 m × 640 m × 1 km, periodic in x and y. Level 0: 20 × 20 × 10 m. Level 1: 10 × 10 × 10 m over x, y = 160–480 m, the full depth |
| Step | 1 s on level 0, two 0.5 s steps on level 1 |
| Air | dry; 300 K up to 300 m, then 10 K/km; 5 m/s geostrophic westerly at 40° N; Smagorinsky LES |
| Sun | 2024-08-05, 15:00–18:00 UTC at 40° N, 100° W (about 08:20–11:20 local solar time) |
| Sky | clear and absorbing: optical depth 0.002 per 10 m layer, no scattering, so no diffuse sky light |
| Ground | two-stream radiation with force-restore balance, as case 2: albedo 0.21, emissivity 0.994, deep temperature 290 K |
| Buildings' balance | prognostic skin, 10 slab layers, convective velocity scale and stability functions on, 16 × 8 rays for the view fractions |

`make_buildings.py` writes the height map `buildings.txt`, the table `buildings.csv` and
`inputs_buildings`, which the deck reads through `FILE =`. Run it with another `--seed` or `--n` for
another set. The refined box holds every building with at least one level-0 cell around it;
`python3 ../../SEB/ibseb_refinement_box.py inputs` confirms it.

**How the faces get their radiation.**
- `erf.ibseb.sun_mode = two_stream`: the faces take the sun's position from the two-stream radiation.
  (Its orbital formula has no equation of time, so the `solar` sun would sit 1.2° away on this date.)
- `erf.ibseb.radiation = two_stream`: each face reads the two-stream column next to it, at its own
  height. It takes the direct beam there (divided by cos z, the direct-normal irradiance), the sky's
  diffuse light and longwave coming down, and the light and longwave coming up from the ground.
- With this sky the beam at height h is S0 cos z exp(-0.002 (100 - h / 10 m) / cos z), with z the
  zenith angle.

## Running it

```
python3 make_buildings.py                # the files are in the directory already; rerun for another set
mpirun -np 4 erf_exec inputs             # about 2 h 30 min on 4 ranks of a laptop
python3 check_random_buildings.py        # exit code = number of failed checks
python3 plot_random_buildings.py         # random_buildings_{map,tskin,theta}.png
```

## What is checked

`check_random_buildings.py` reads `buildings.csv`, the building report `ibseb_buildings.csv`, the
face dumps (`faces/set.*` on level 0, `faces/set.lev1.*` on level 1), `radiation_diag.csv` and the
last 2D plotfile. Each check, and what it would catch:

1. **Buildings.** Both levels hold all 8 buildings, numbered alike, their roofs centred on
   `buildings.csv`. Catches a level that loses or renumbers a building.
2. **Energy balance.** Every building closes its balance on both levels (residual below 1e-3 W/m²).
3. **The sun and the beam.** On every sunlit, unshaded roof of level 0, at every face dump, the
   top-of-atmosphere irradiance implied by the roof's direct sunlight at its own height equals the
   sweep's, to 1e-5. There is no diffuse sky, and the sun is east of the meridian. Catches a face
   reading another height or another sun. (A dump carries the time at the end of its step and the
   sun of its start, so it is compared with the sweep one step earlier.)
4. **The sun climbs.** The shaded share of the face area falls on both levels, from the first report
   with radiation (the report at the start comes before any sweep).
5. **Warming.** Every building warms on both levels.
6. **Levels agree.** Each building's top roof within 0.5 K and its walls within 1 K between the
   levels, at the last dump. The whole-building means are reported, not checked (see below).
7. **Open ground.** The ground away from the footprints warms on both levels.

The checker also prints each building's top roof, ledges and walls on both levels: area, absorbed
sunlight and skin temperature. The results below rest on those numbers.

## Results

These numbers come from the run before the last change to the faces (each wall now reads the
light and longwave coming up from the ground at its own height, not at the ground). That change
moves the walls' incoming radiation by up to about 4 W/m² on the 60 m walls; a rerun will update
the tables below.

On 4 MPI ranks of a laptop shared with other runs, the 10,800 level-0 steps took about 2 h 30 min.
The building balance cost 0.25 ms per step on level 0 and 0.11 ms on level 1, on the slowest rank.

All 7 checks pass:

| Check | Result |
|---|---|
| 1. Buildings | 8 and 8, roof centres on `buildings.csv` |
| 2. Energy balance | largest residual 3.3e-8 W/m² over 2,896 rows |
| 3. The sun and the beam | 11,034 sunlit roof faces at 5 heights, largest difference 4.8e-7; no diffuse sky; azimuth 99.5° to 156.6° |
| 4. The sun climbs | zenith 52.7° to 24.8°; shaded face area 5.3 % to 3.2 % on level 0, 6.5 % to 4.2 % on level 1 |
| 5. Warming | final means 314.0–321.7 K on level 0, 312.6–319.6 K on level 1 (from 300 K) |
| 6. Levels agree | level 1 minus level 0: top roofs -0.12 to -0.01 K, walls -0.64 to +0.33 K |
| 7. Open ground | 310.52 K in the refined area (level 1), 310.42 K outside it (level 0) at 18:00 UTC |

The sun the faces saw, with the mean direct-normal irradiance over the level-0 faces at their own
heights:

| UTC | Zenith | Azimuth | Direct normal |
|---|---|---|---|
| 15:01 | 52.5° | 99.7° | 958 W/m² |
| 16:00 | 41.6° | 112.1° | 1017 W/m² |
| 17:00 | 31.8° | 129.8° | 1050 W/m² |
| 18:00 | 24.8° | 156.6° | 1065 W/m² |

Each building's parts at 18:00 UTC, level 0 / level 1: the absorbed sunlight and the skin
temperature of its top roof and of its walls.

| Building | Height, material | Top roof: SW [W/m²], skin [K] | Walls: SW [W/m²], skin [K] | Mean skin [K] |
|---|---|---|---|---|
| 1 | 60 m, timber | 586 / 586, 324.0 / 324.0 | 122 / 124, 316.5 / 315.9 | 318.1 / 316.4 |
| 2 | 30 m, timber | 582 / 582, 324.7 / 324.7 | 137 / 142, 318.9 / 319.0 | 321.7 / 319.6 |
| 3 | 60 m, brick | 683 / 683, 323.1 / 323.0 | 149 / 149, 312.0 / 311.6 | 314.0 / 312.6 |
| 4 | 20 m, brick | 677 / 677, 323.3 / 323.3 | 163 / 170, 311.5 / 311.8 | 317.3 / 316.8 |
| 5 | 10 m, concrete | 724 / 724, 321.8 / 321.7 | 172 / 172, 310.5 / 310.8 | 317.7 / 317.5 |
| 6 | 50 m, timber | 584 / 584, 324.1 / 324.0 | 142 / 146, 317.9 / 317.5 | 318.3 / 318.1 |
| 7 | 50 m, timber | 584 / 584, 325.5 / 325.5 | 132 / 133, 317.2 / 316.8 | 318.5 / 317.9 |
| 8 | 20 m, brick | 677 / 677, 323.4 / 323.3 | 155 / 156, 312.3 / 312.6 | 316.8 / 316.3 |

What the numbers show:

- **Roofs absorb four to five times what walls do.** The top roofs absorb 582–724 W/m² at 18:00
  UTC and the walls 122–172 W/m². Buildings that are mostly roof absorb the most per square metre:
  the 10 m concrete block (5) about 525 W/m² over all its faces, the 60 m towers (1, 3) 182–243.
- **A roof's sunlight depends on its material and a little on its height.** Timber roofs absorb
  582–586 W/m², brick 677–683 and concrete 724. Higher roofs absorb slightly more, as less
  absorbing air lies above them.
- **The two levels agree on roofs and walls.** Their skin temperatures agree to 0.12 K on the top
  roofs and 0.64 K on the walls.
- **They differ at the buildings' edges.** Each level draws a building's edge as a step one cell
  wide, which faces upward (a ledge). Level 0 has one ledge 20 m wide, level 1 narrower ones at
  several heights. On buildings 1–3 the level-1 ledges absorb 327–395 W/m² against 483–579 on
  level 0, and run 4–9 K cooler. That makes each building's mean skin 0.2–2.2 K cooler on level 1.
- **The open ground is the same on both levels:** 310.52 K in the refined area, 310.42 K outside it.
  Under and beside the footprints the ground runs hotter (312.09 K on level 1), as it still takes
  the sun there.

## Compared with the clear-sky faces

The run before this option, with the clear-sky formulas set to match the same beam (a gray sky of
emissivity 0.83, a fixed 300 K ground), is the reference here. Both runs are from before the last
change to the faces.

- **Sunlight agrees.** The absorbed sunlight agrees within 1 % on every building (largest 0.92 %).
- **The buildings end warmer.** On level 0 the buildings end 1.2–4.2 K warmer and their walls
  1.8–5.2 K warmer; on level 1, 1.2–5.0 K and 1.9–5.8 K.
- **Why: more longwave from the sky and the ground.** At 18:00 UTC the walls take in about 47 W/m²
  more longwave than with the clear-sky formulas. About 19 W/m² comes from the sky: the column's sky
  is warmer than the gray sky (an effective emissivity near 0.93, from case 2's longwave optical
  depth). About 29 W/m² comes from the ground: the walls see open, sunlit ground near 312.7 K,
  where the clear-sky run had 300 K. The roofs see no ground and take about 40 W/m² more, all from
  the sky.

## Plots

`plot_random_buildings.py` draws:

- **`random_buildings_map.png`:** the last 2D plotfile from above: the open ground's skin, the
  footprints and the level-1 roofs.
- **`random_buildings_tskin.png`:** the skin temperature of every building on both levels against
  time, with the open ground.
- **`random_buildings_theta.png`:** the potential temperature 15 m up on level 1, with the wind.

## Limits and notes

- **The ground under and beside the buildings.** The two-stream columns pass through the buildings,
  so the ground balance still covers the footprints and takes the sun there. The checks leave the
  footprints out. The ground beside a wall is open, sunlit ground, which the walls see.
- **Sunlight counted twice.** Where the buildings stand, the faces and the ground both absorb the
  sunlight: about 7 % more than the domain receives, and 22 % more within the refined
  area, at 18:00 UTC.
- **Mixed-layer depth.** The convective velocity scale (`erf.ibseb.z_i_mode = bulk_ri`) takes the
  mixed-layer depth from level 0's profile on both levels.
- **Ground temperature output.** The two-stream diagnostics file leaves `T_s_mean` empty (`nan`) with
  the prognostic balance; the ground temperatures above come from the 2D plotfiles.
