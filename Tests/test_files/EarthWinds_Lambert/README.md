# Earth-relative winds on a Lambert-conformal grid

## Purpose

Checks the earth-relative winds `u_east`, `v_north`, `wind_speed` and `wind_dir` on a
grid that is rotated with respect to east and north, both as station columns and as 2D
sampled-level fields. The 25 km Lambert-conformal grid (TRUELAT1/2 = 30/60,
STAND_LON = -120) turns by up to 11.4 degrees from east at its corners, so a field that is
not rotated, or is rotated the wrong way, is off by up to 0.2 rad of the wind.

The reference is the closed-form Lambert rotation,

    alpha = -n (lon - STAND_LON),   n = ln(cos T1 / cos T2) / ln(tan(45 - T1/2) / tan(45 - T2/2)),

which matches WRF's own `COSALPHA`/`SINALPHA` in this file to 4e-6 degrees. ERF measures
the angle from its lat/lon arrays instead, so the comparison is independent.

## Files

- `EarthWinds_Lambert.i`: one step from `wrfinput_d01` (no lateral boundary file:
  outflow sides). Writes three stations (NW and NE placed by `lat`/`long`, SW by `x`/`y`)
  and a 2D plotfile with `lon_m` and the sampled `u_east`/`v_north` at `k = 10`.
- `Tests/check_station_earth_winds.py`: per station, the header angle against the
  expected one (0.01 degrees), and the rotation, speed and direction of every row.
- `Tests/check_2d_earth_winds.py`: every cell of the 2D fields against the 3D
  `x_velocity`/`y_velocity` rotated by the closed-form angle (0.003 rad of the wind).

## Data and acknowledgements

The case needs one file that is not in the repository: `wrfinput_d01`, a WRF real.exe
(V4.8.0) initial-condition file for 2024-08-05 00 UTC on a 100 x 100 x 79 grid at 25 km
(Lambert conformal, `CEN_LAT` = 33, `CEN_LON` = -120, MODIS land use). It is about 24 MB
and is not redistributed here. Place it at

    ${ERF_TEST_EXTRA_FILES_DIRECTORY}/EarthWinds_Lambert/wrfinput_d01

Any WRF file on a Lambert grid can be used in its place, with the stations' expected
angles and the projection constants in `Tests/CTestList.cmake` changed to match.

## How to run

Configure with `-DERF_ENABLE_NETCDF=ON` (parallel NetCDF) and, if the file is not under
`Tests/test_files/`, `-DERF_TEST_EXTRA_FILES_DIRECTORY=<dir>`. Then

    ctest -R EarthWinds_Lambert

Without the file, CMake prints `Skipping EarthWinds_Lambert` and adds no test.

## Expected Results

| Check | Value |
|---|---|
| NW / NE / SW grid angle (ERF vs closed form) | 8.7840 / 8.7845, -7.6401 / -7.6411, 6.7328 / 6.7324 degrees |
| Station rotation, speed, direction | relative error below 1e-9 |
| 2D `u_east`/`v_north`, worst cell | 2.0e-3 rad of the wind (edge columns, one-sided difference); interior about 3e-5 |
| 2D fields on `development` before #4700 | 0.20 rad (not rotated): the check fails |

## CTest

`EarthWinds_Lambert` (labels `regression`, `station`), added only when NetCDF is enabled
and the file above exists.
