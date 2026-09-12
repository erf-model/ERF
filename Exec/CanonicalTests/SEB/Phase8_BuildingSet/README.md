# SEB/Phase8_BuildingSet

Canonical case for the balance on a set of buildings from a height map:
mutual shadowing, the building part of the view fractions, several
materials, the per-building report, and the cost of the face list, with
the immersed forcing's wall law snapped to whole cells and the wall
function beyond neutral switched on.

```
./run_buildingset.sh /path/to/erf_exec                   # NP=4, about an hour
python3 plot_buildingset.py ibseb_set.csv faces/set --plotfile plt43200
```

## The scenario

Four buildings on a 480 m periodic domain at 10 m, from a nodal height map
(`make_buildings.py`), numbered as the balance numbers them (scan order):

1. a 60 m concrete slab, 30 x 60 m, at x = 150-180, y = 200-260;
2. a 20 m timber block, 40 x 40 m, north of the slab at y = 300-340;
3. a 40 m brick cube, 40 x 40 m, east of the slab at x = 230-270, in
   its morning shadow;
4. a 20 m timber block off on its own at x = 320-360, y = 150-190.

Three materials from `materials.csv` by building; a 3 m/s westerly,
neutral at 300 K with a MOST ground; Boulder on the June solstice from
05:00 solar time for six hours with the prescribed clear-sky provider and
a gray sky; the convective velocity scale (bulk Richardson depth) and the
stability functions on; `erf.if_snap_partial_cells = true`, which a
height-map set needs (`Exec/RegTests/ImmersedForcingTest/PartialCells`).

## What is checked

1. Four buildings with the expected face counts and materials.
2. The balance residual stays below 1e-3 W/m2 on every building.
3. Mutual shadowing: at sunrise the slab is the most shaded building,
   because the 40 m cube east of it throws its shadow onto the slab's
   east wall (the slab shades the cube in the afternoon, outside this
   run), and it clears by late morning; the 20 m blocks are free of
   shadow by then. Every building also shades its own stepped rim at low
   sun.
4. The building view fraction: the facing walls of the slab and the cube
   see more building than the far block's walls (every wall sees some,
   its own stepped rim and, across the periodic boundary, the others).
5. Materials show: the timber roofs end warmer than the concrete roof.
6. The wall function beyond neutral is active: w* positive on the sunlit
   faces, most roofs unstable by the end (the shaded rim roofs may be
   slightly stable).
7. Timing: faces per rank and the cost of the balance per step, reported
   for estimating a city-scale case.

## What happens through the morning

- **Sunrise (04:35).** The run starts at 05:00 with the sun 9 degrees up in
  the east-north-east. The cube's shadow lies across the slab's east wall,
  so the slab is the most shaded building (17 percent of its faces, 6
  percent by late morning); the 20 m blocks shade only their own stepped
  rims, and are free of shadow after 10:30. The cube keeps a little
  self-shadow at 11:00 from its 20 m rim step.
- **Materials.** The two timber blocks (light cladding, 15 cm) warm
  fastest and end at 332.7 and 332.5 K on their roofs; the brick cube
  (25 cm) reaches 324 K and the concrete slab (30 cm) 320.5 K. Their
  mean skin temperatures rise together at first and split by material
  after 06:00, the identical timber blocks tracking each other to within
  0.2 K although one sits north of the slab and one on its own.
- **View fractions.** The slab's east wall sees 28 percent building, the
  cube's west wall 48 percent (the taller slab fills its view), the far
  block's walls 18 percent (their own rim and the periodic neighbours).
- **The wall function.** The convective scale is 0.1 to 0.7 m/s on the 268
  sunlit faces at 11:00, and 93 percent of the roofs are unstable with
  their own Obukhov length; the shaded rim roofs run slightly stable.
- **Cost.** 0.4 ms per step on the slowest rank for 157 faces, 6 ms to
  build the list and sample the view fractions; the balance is a
  negligible part of the 126 minutes the six hours took on four ranks.
- **The atmosphere.** The height-map set runs with the immersed forcing
  snapped to whole cells and the implicit drag; without the snap the slab
  drives the density negative at two minutes (see the plan file's
  findings), and the balance's residual stays at 3e-8 W/m2 throughout.

## Reference output (4 ranks)

```
== building set, 6 h from 05:00 (4 ranks)
  four buildings with their materials: PASS (building 1: 172 faces, material 1, building 2: 72 faces, material 3, building 3: 116 faces, material 2, building 4: 72 faces, material 3)
  balance residual on every building all morning: PASS (max 3.4e-08 W/m2 over 1444 rows)
  at sunrise the slab is the most shaded building (the cube's shadow on its east wall) and clears by late morning: PASS (slab 0.168 at 05:30-06:10 vs 0.061 after 10:30; others early 0.108, 0.084, 0.112)
  the 20 m blocks are free of shadow by late morning: PASS (north block 0.000, far block 0.000 after 10:30)
  the facing walls of the slab and the cube see more building than the far block's walls: PASS (slab east 0.28, cube west 0.48, far block walls 0.179)
  timber roofs end warmer than the concrete roof (light cladding warms faster): PASS (roof means at 11.0 h: slab 320.5, north block 332.7, cube 324.2, far block 332.5 K)
  w* positive on the sunlit faces and most roofs unstable by the end (the shaded rim roofs may be stable): PASS (w* 0.08-0.72 m/s on 268 sunlit faces, 93 % of 140 roofs unstable)
  cost line reported: PASS (lev=0 advance_ms_per_step_max=0.3804 faces_per_rank_max=157 ranks=4 init_s=0.005792)
building set: PASS
ALL PASS
```
