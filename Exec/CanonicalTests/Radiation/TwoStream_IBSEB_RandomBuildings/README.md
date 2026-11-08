# Two-stream radiation over random buildings, on two levels

The third of the two-stream test cases. It takes the ground of case 2
(`../TwoStream_NoahMP_vs_ForceRestore`, the two-stream radiation with the force-restore surface
balance) and stands a seeded random set of buildings on it. The two-stream radiation heats the air
and the ground. The immersed-boundary surface energy balance (`erf.ibseb`) heats the walls and roofs,
and its faces take their radiation from the same two-stream columns. The columns see the buildings
(`erf.ibseb.two_stream_buildings`, on by default): no radiative heating inside them, the ground
under them held, and the ground they shade without the direct sunlight. The buildings sit on a
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
- With this sky the beam at height h is S0 cos z exp(-0.002 (100 - h / 10 m) / cos z), with S0 the
  sun's irradiance at the top of the atmosphere and z the zenith angle.

## Running it

```
python3 make_buildings.py                # the files are in the directory already; rerun for another set
mpirun -np 4 erf_exec inputs             # about 2 h 20 min on 4 ranks of a laptop
python3 check_random_buildings.py        # exit code = number of failed checks
python3 plot_random_buildings.py         # random_buildings_{map,tskin,theta}.png
```

## What is checked

`check_random_buildings.py` reads `buildings.csv`, the building report `ibseb_buildings.csv`, the
face dumps (`faces/set.*` on level 0, `faces/set.lev1.*` on level 1), `radiation_diag.csv` and the
first and last 2D plotfiles. Each check, and what it would catch:

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
8. **Ground under the buildings.** Every column under a footprint ends at the skin temperature it
   started with, on both levels, and absorbs no sunlight. Catches a footprint that still takes the
   sun (in the run before this check, by up to 16 K in 40 minutes). The checker also reports the
   shaded share of the open ground and the sunlight it absorbs.

The checker also prints each building's top roof, ledges and walls on both levels: area, absorbed
sunlight and skin temperature. The results below rest on those numbers.

## Results

Two runs of the deck on development ae7ef9d47 plus this change, on 4 MPI ranks of a laptop shared
with other runs: the option on (the default) and off (`erf.ibseb.two_stream_buildings = false`, the
behaviour of #4232). Each took about 2 h 20 min for the 10,800 level-0 steps. The building balance
cost 0.20 ms per step on level 0 and 0.11 ms on level 1, on the slowest rank. (The numbers differ
from those in the #4232 run of this deck: the dynamics changed between the two, and the deck now
pins the acoustic substep count, see `inputs`.)

With the option on all 8 checks pass. With it off, checks 1 to 7 pass and check 8 fails, as it
should: the footprints absorb the full 761 W/m² and no open ground is shaded.

| Check | On (default) | Off |
|---|---|---|
| 1. Buildings | 8 and 8, roof centres on `buildings.csv` | same |
| 2. Energy balance | largest residual 3.3e-8 W/m² over 2,896 rows | same |
| 3. The sun and the beam | 11,034 sunlit roof samples at 5 heights, largest difference 4.8e-7; no diffuse sky; azimuth 99.5° to 156.6° | same |
| 4. The sun climbs | zenith 52.7° to 24.8°; shaded face area 5.3 % to 3.2 % on level 0, 6.5 % to 4.2 % on level 1 | same |
| 5. Warming | final means 313.0–319.9 K on level 0, 311.5–316.9 K on level 1 (from 300 K) | 313.3–321.2 K, 311.9–319.0 K |
| 6. Levels agree | level 1 minus level 0: top roofs -0.13 to +0.10 K, walls -0.98 to -0.03 K | -0.06 to +0.08 K, -0.54 to +0.28 K |
| 7. Open ground | 309.67 K in the refined area (level 1), 309.66 K outside it (level 0) at 18:00 UTC | 309.76 K, 309.65 K |
| 8. Ground under the buildings | 95 footprint columns on both levels at their starting skin, 0 W/m² absorbed; 3.7 % of level 1's open ground shaded, 0 W/m² there against 761 W/m² in the sun | **fails:** footprints absorb 761 W/m², 0 % shaded |

The sun the faces saw, with the mean direct-normal irradiance over the level-0 faces at their own
heights (the same in both runs):

| UTC | Zenith | Azimuth | Direct normal |
|---|---|---|---|
| 15:01 | 52.5° | 99.7° | 958 W/m² |
| 16:00 | 41.6° | 112.1° | 1017 W/m² |
| 17:00 | 31.8° | 129.8° | 1050 W/m² |
| 18:00 | 24.8° | 156.6° | 1065 W/m² |

The buildings at 18:00 UTC with the option on, level 0 / level 1: mean skin temperature, absorbed
sunlight (SW), net longwave (LW) and sensible heat flux (H) over all faces. The absorbed sunlight is
the same with the option off (the faces' light comes from the columns at their own heights either
way); the other three columns change because the walls see a different ground.

| Building | Height, material | Mean skin [K] | SW [W/m²] | Net LW [W/m²] | H [W/m²] |
|---|---|---|---|---|---|
| 1 | 60 m, timber | 316.4 / 314.3 | 225 / 181 | -84 / -75 | 105 / 82 |
| 2 | 30 m, timber | 319.9 / 316.9 | 306 / 237 | -94 / -81 | 174 / 129 |
| 3 | 60 m, brick | 313.0 / 311.5 | 243 / 206 | -61 / -54 | 71 / 55 |
| 4 | 20 m, brick | 317.2 / 316.4 | 408 / 378 | -90 / -82 | 158 / 139 |
| 5 | 10 m, concrete | 317.3 / 316.9 | 525 / 509 | -99 / -95 | 183 / 178 |
| 6 | 50 m, timber | 317.2 / 316.6 | 230 / 221 | -85 / -85 | 102 / 91 |
| 7 | 50 m, timber | 316.9 / 315.8 | 219 / 204 | -89 / -83 | 94 / 84 |
| 8 | 20 m, brick | 316.2 / 315.4 | 356 / 330 | -81 / -74 | 135 / 121 |

The same buildings with the option off:

| Building | Mean skin [K] | SW [W/m²] | Net LW [W/m²] | H [W/m²] |
|---|---|---|---|---|
| 1 | 317.2 / 315.5 | 225 / 181 | -75 / -64 | 116 / 94 |
| 2 | 321.2 / 319.0 | 306 / 237 | -96 / -80 | 159 / 110 |
| 3 | 313.3 / 311.9 | 243 / 206 | -52 / -41 | 93 / 80 |
| 4 | 316.8 / 316.3 | 408 / 378 | -87 / -75 | 180 / 161 |
| 5 | 317.0 / 316.7 | 525 / 509 | -96 / -91 | 225 / 221 |
| 6 | 317.8 / 317.6 | 230 / 221 | -78 / -75 | 115 / 106 |
| 7 | 317.0 / 316.3 | 219 / 204 | -78 / -69 | 114 / 106 |
| 8 | 316.5 / 316.0 | 356 / 330 | -79 / -67 | 137 / 121 |

Each building's top roof and walls at 18:00 UTC with the option on, level 0 / level 1: absorbed
sunlight and skin temperature.

| Building | Top roof: SW [W/m²], skin [K] | Walls: SW [W/m²], skin [K] |
|---|---|---|
| 1 | 586 / 586, 324.6 / 324.6 | 121 / 123, 314.4 / 313.4 |
| 2 | 582 / 582, 322.9 / 323.0 | 137 / 142, 316.7 / 315.8 |
| 3 | 683 / 683, 322.2 / 322.1 | 149 / 148, 310.9 / 310.3 |
| 4 | 677 / 677, 322.7 / 322.7 | 163 / 169, 311.5 / 311.4 |
| 5 | 724 / 724, 321.2 / 321.2 | 172 / 172, 310.2 / 310.2 |
| 6 | 584 / 584, 324.4 / 324.3 | 142 / 146, 316.4 / 315.5 |
| 7 | 584 / 584, 324.7 / 324.7 | 131 / 133, 315.3 / 314.4 |
| 8 | 677 / 677, 322.7 / 322.7 | 155 / 156, 311.8 / 311.4 |

What the numbers show:

- **What the option changes is the ground, and through it the walls.** With it off the ground under
  and beside the footprints takes the sun and ends at 311.5 K on level 1, 1.7 K above the open
  ground; with it on it is held at 306.7 K (the footprints at their 300 K start, their shaded
  surroundings cooler than the open ground). The open ground away from the buildings is the same
  within 0.1 K either way, and so is the sunlight every face absorbs.
- **The walls end cooler with the option on.** On level 0 the buildings' mean skin is 0.3 to 1.3 K
  lower (0.2 to 2.1 K on level 1): the walls see cooler ground, and their net longwave is 2 to
  12 W/m² more negative. The sensible heat flux of the low buildings (4, 5, 8) drops by 22 to
  42 W/m²; the air near them is no longer heated by sunlit ground under and around the footprints.
  Building 2's heat flux rises by 15 W/m² while its skin falls, so the air in its street must be
  cooler with the option on; this case does not separate the two effects further.
- **Roofs absorb four to five times what walls do.** The top roofs absorb 582–724 W/m² at 18:00
  UTC and the walls 121–172 W/m². Buildings that are mostly roof absorb the most per square metre:
  the 10 m concrete block (5) 525 W/m² over all its faces, the 60 m towers (1, 3) 181–243.
- **A roof's sunlight depends on its material and a little on its height.** Timber roofs absorb
  582–586 W/m², brick 677–683 and concrete 724. Higher roofs absorb slightly more, as less
  absorbing air lies above them.
- **The two levels agree on roofs and walls.** Their skin temperatures agree to 0.13 K on the top
  roofs and 0.98 K on the walls with the option on.
- **They differ at the buildings' edges.** Each level draws a building's edge as a step one cell
  wide, which faces upward (a ledge). Level 0 has one ledge 20 m wide, level 1 narrower ones at
  several heights. On buildings 1–3 the level-1 ledges absorb 327–395 W/m² against 483–579 on
  level 0, and run 4–9 K cooler. That makes each building's mean skin 0.4–3.0 K cooler on level 1.
- **The shadow on the ground is small at these sun heights:** 3.7 % of level 1's open ground at
  18:00 UTC with the sun 65° up, in the shaded columns the ground absorbs no sunlight against
  761 W/m² in the sun.

## Compared with the clear-sky faces

The numbers in this section are from the #4232 run of this deck (the option did not exist yet; the
faces saw the same sunlight and the ground under the buildings took the sun).

The run before this option, with the clear-sky formulas set to match the same beam (a gray sky of
emissivity 0.83, a fixed 300 K ground), is the reference here.

- **Sunlight agrees.** The absorbed sunlight agrees within 1 % on every building (largest 0.81 %).
- **The buildings end warmer.** On level 0 the buildings end 1.2–3.9 K warmer and their walls
  1.8–4.9 K warmer; on level 1, 1.3–4.5 K and 1.8–5.3 K.
- **Why: more longwave from the sky and the ground.** At 18:00 UTC the level-0 walls take in about
  46 W/m² more longwave than with the clear-sky formulas. About 19 W/m² of it comes from the sky,
  which is warmer in the two-stream column than the gray sky (an effective emissivity near 0.93,
  from case 2's longwave optical depth). About 28 W/m² comes from the ground: the walls see open,
  sunlit ground near 312 K, where the clear-sky run had 300 K. The roofs see no ground and take
  about 40 W/m² more, all from the sky.

## Plots

`plot_random_buildings.py` draws:

- **`random_buildings_map.png`:** the last 2D plotfile from above: the open ground's skin, the
  footprints and the level-1 roofs.
- **`random_buildings_tskin.png`:** the skin temperature of every building on both levels against
  time, with the open ground.
- **`random_buildings_theta.png`:** the potential temperature 15 m up on level 1, with the wind.

## Limits and notes

- **The ground under and beside the buildings.** The ground under a building is held and the
  shaded ground gets no direct sunlight, but the two-stream columns still pass through the
  buildings: the light the column's ground reflects, which the walls see, still includes the beam
  the shaded ground no longer absorbs, and the ground beside a wall still sees the whole sky.
- **Mixed-layer depth.** The convective velocity scale (`erf.ibseb.z_i_mode = bulk_ri`) takes the
  mixed-layer depth from level 0's profile on both levels.
- **Ground temperature output.** The two-stream diagnostics file leaves `T_s_mean` empty (`nan`) with
  the prognostic balance; the ground temperatures above come from the 2D plotfiles.
