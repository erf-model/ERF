# Flood engine example

## Purpose

The smallest complete flood engine directory, built from outside ERF's `Source` tree.
It checks the flood coupling plumbing (`Source/FloodCoupling`): the out-of-tree build
(CMake `ERF_FLOOD_ENGINE_DIRS`, GNU make `FLOOD_HOME`), the engine registry, the
precipitation handed to the engine each step, the flood plotfile, the flood fields in
the ERF plotfile and checkpoint/restart.

The engine is a bucket: each flood cell keeps all the water that falls on it, with no
runoff, infiltration or evaporation,

    h += (P_liquid + P_frozen) / rho_w,    rho_w = 1000 kg/m3,

with the step's precipitation `P` in kg/m2.

## Files

| File | What it is |
|------|------------|
| `ERF_FloodEngineExample.cpp` | The engine, registered as `example` |
| `ERFFloodEngine.cmake` | CMake hook (`erf_flood_engine_add`, `erf_flood_engine_add_unit_tests`) |
| `Make.package` | GNU make hook |
| `ERF_GTestFloodEngineExample.cpp` | Unit test: the engine is registered and declares `flood_depth` |
| `inputs_flood_example` | Kessler rain in a laterally periodic layer, 60 s, with a cool bubble off the centre so the rain varies about threefold over the columns |
| `check_depth.py` | Compares every flood cell with the rain of its own column; uses `Exec/Tools/erf_plotfile.py` |
| `sounding_moist` | The moist sounding (a copy of `Tests/test_files/RegridKeepsPrecipAccum/sounding_moist`) |

## Data

None. The sounding is an idealized profile from ERF's own test files.

## How to run

```bash
cmake -S . -B build-flood -DERF_ENABLE_MPI=ON -DERF_ENABLE_TESTS=ON \
      -DERF_ENABLE_FLOOD=ON -DERF_FLOOD_ENGINE_DIRS=$PWD/Tests/FloodEngineExample
cmake --build build-flood -j 8
cd build-flood && ctest -R Flood_ExampleEngine --output-on-failure
```

With GNU make: `make USE_FLOOD=TRUE FLOOD_HOME=<ERF>/Tests/FloodEngineExample` in an `Exec` case directory.

## Expected results

After 120 steps of 0.5 s the surface rain ranges from about 4.9e-4 to 1.5e-3 kg/m2. On every
flood cell, `flood_depth` times 1000 kg/m3 equals `precip_rain_accum` of the column above it to
round-off (worst relative error about 1e-15; `check_depth.py` passes below 1e-9 and fails when
the rain varies less than 1.5 times, which would hide a wrong column). A restart from
`chk00060` reproduces the straight run's plotfiles exactly.

## CTest

All need 2 ranks, double precision and `fcompare`.

| Test | What it checks |
|------|----------------|
| `Flood_ExampleEngine` | The straight run, the per-cell depth check, `flood_depth` in the ERF plotfile, the restart comparison, and a restart onto another `erf.flood.grid_ratio` that must stop |
| `Flood_ExampleEngine_SpinUp` | A flood restarted from an atmosphere-only checkpoint at t = 30 s holds the rain of the second half only |
| `Flood_ExampleEngine_FinerLevel` | The flood grid on a static level-1 box over 400 to 1200 m (an offset region) |
| `Flood_ExampleEngine_UnderFinerLevel` | The flood on level 0 under that box: a restart reproduces both plotfiles, with the correction ERF averages down from level 1 kept in the checkpoint |
| `Flood_Abort_*` | The start-up refusals, each naming the input to change |
