# GeoWindTable_TwoLevel

## Purpose

Checks the height-varying geostrophic wind read from `erf.abl_geo_wind_table`
(erf-model/ERF#4605, #4606): on a refined level that is also refined in z, and on
tables that cannot be used.

## Files

| file | role |
| --- | --- |
| `GeoWindTable_TwoLevel.i` | 800 x 800 x 1000 m, 8 x 8 x 16 cells, at rest, Coriolis factor 1e-4 1/s at the pole; level 1 covers the whole horizontal extent and is refined by 2 in x, y and z; one 1 s step |
| `geo_wind_table` | Ug = 2 + z/100 m/s, Vg = 0, two rows and a comment line |
| `geo_wind_table_empty` | a comment and a blank line, no rows |
| `geo_wind_table_short_row` | a second row with two numbers (rows must hold exactly three) |
| `geo_wind_table_decreasing` | heights 0, 600, 500, 1000 m |
| `../../check_geowind_profile.py` | the check for the two-level run (Python standard library only) |

## How to run

```bash
mpiexec -n 2 erf_exec GeoWindTable_TwoLevel.i
python3 ../../check_geowind_profile.py plt00001
```

## Expected results

Over the first step the only force on the fluid at rest is the geostrophic one, so
`v = f dt Ug(z)` at every cell centre; the neglected terms are of relative size
`(f dt)^2 = 1e-8`. The check passes below 1e-6 on both levels (1e-4 in a single-precision build); the
code gives 1.7e-9.
Before the fix, level 1 read level 0's z levels, past their end on a level refined in
z, and the run stopped with a signal in the first step.

The three other tables stop at start-up with a message naming the file and the line.
Before the fix the empty table read past the end of an empty list, the short row
reused the previous row's Vg, and the decreasing heights were interpolated as given.

A fine level that does not cover the whole horizontal extent adds its own small
velocity changes at its side faces in the first step (5 to 13 % of v here), so the
deck refines the whole horizontal extent.

## CTest

`ctest -R GeoWindTable` runs `GeoWindTable_TwoLevel` (needs Python 3) and the three
abort tests `GeoWindTable_Empty`, `GeoWindTable_ShortRow` and `GeoWindTable_Decreasing`.
