# Two-stream radiation with Noah-MP versus the force-restore balance, on two levels

The two surface models the two-stream radiation can drive, run side by side on the same
refined grid, the same atmosphere and the same sun, so that the simpler one (the
force-restore surface energy balance) can be judged against the fuller one (Noah-MP).

## Setup

Shared by both cases (`inputs_common`):

- **Grid.** 2 km x 2 km x 4 km, 8 x 8 x 64 cells (62.5 m in z). Level 1 covers the middle 1 km x 1 km,
  refined by 2 in x and y, spanning z, so it runs the two-stream sweep on its own columns.
  It stays fixed for the run.
- **Time step.** `erf.cfl = 0.5` on every level.
- **Atmosphere.** A 300 m mixed layer at 300 K under a 6 K/km inversion, 14 g/kg at the
  surface drying to 8 g/kg at 300 m and 2 g/kg at 4 km, a 5 m/s geostrophic westerly at 40 N, MRF boundary
  layer, Kessler microphysics without rain.
- **Sun.** 2024-08-05 from 15:00 to 21:00 UTC at 40 N, 100 W (about 08:20 to 14:20 local
  solar time).

The two cases differ only in the surface:

| | `inputs_noahmp` | `inputs_force_restore` |
|---|---|---|
| Surface model | Noah-MP on both levels | force-restore balance on both levels |
| Land | grassland (IVGTYP 10) on silty clay loam, 290 K soil at 0.25 m3/m3; level 1 has its own land setup file | none: the balance's own skin |
| Surface albedo and emissivity | Noah-MP's: emissivity 0.994; albedo 0.229 at 15:10 UTC falling to 0.201 at 18:00 | 0.21 (about Noah-MP's daytime mean) and 0.994 |
| Surface humidity | Noah-MP's evapotranspiration | fixed at the sounding's 14 g/kg |
| Deep temperature | Noah-MP's soil column | restore to 290 K (the soil's) over a day |
| Land step | 60 s (`NOAH_TIMESTEP`) | every ERF step |
| H and LE | Noah-MP's, applied through the surface layer | the surface layer's, which the balance removes (`seb_turbulent_flux_source = surface_layer`), with the surface layer taking its surface temperature from the balance's skin (`seb_surface_layer_uses_skin`) |

Noah-MP runs on both levels because `namelist.erf` names a land setup file for each
(`ERF_SETUP_FILE_01` and `_02`). `make_land_files.py` writes them; level 1's places itself
in level 0 with the WRF attributes `I_PARENT_START = J_PARENT_START = 3` and
`PARENT_GRID_RATIO = 2`.

## Running it

Needs an executable built with `-DERF_ENABLE_NOAHMP=ON` (and so a parallel NetCDF), `ncgen`
on the `PATH`, and Python 3 with numpy, matplotlib and yt:

```
./run_comparison.sh /path/to/erf_exec [nranks]
```

It writes the land files, runs both cases under `runs/`, and calls `compare.py`, which
writes `comparison.csv` (level means every 10 minutes) and `comparison.png`.

## Results

From the release build on 2 MPI ranks (about 11 minutes per case). Level means, Noah-MP
vs force-restore:

| | 18:50 UTC (3.8 h, Noah-MP's peak H) | 21:00 UTC (6 h) |
|---|---|---|
| Absorbed shortwave [W/m²] | 780.6 vs 771.6 | 633.7 vs 636.8 |
| Net longwave, down [W/m²] | -123.2 vs -76.3 | -118.3 vs -75.9 |
| Skin temperature [K] | 313.4 vs 304.9 | 314.0 vs 305.5 |
| H [W/m²] | 282.5 vs 95.7 | 225.5 vs 61.8 |
| LE [W/m²] | 166.5 vs 321.8 | 158.7 vs 293.2 |
| PBL height [m] | 1285 vs 863 | 1595 vs 1034 |

- **The radiation agrees.** The shortwave each surface absorbs differs by about 1 %; the
  remaining difference is the albedo, which Noah-MP varies with the sun and the balance holds
  at 0.21.
- **The partition does not.** The force-restore surface puts about twice Noah-MP's energy into
  evaporation and a third of its into H, so its skin stays 8.5 K cooler, loses 45 W/m² less
  longwave, and grows a boundary layer two-thirds as deep. The cause is the surface humidity.
  The balance has no soil water: the surface layer evaporates from a fixed surface mixing
  ratio, and ERF requires that one to equal the sounding's surface value (14 g/kg), so the
  surface behaves as if wet. Noah-MP limits evapotranspiration through its soil moisture,
  stomata and canopy. The balance's own prognostic surface humidity (`seb_q_sfc`) is not
  handed to the surface layer: `seb_surface_layer_uses_skin` couples the temperature only.
- **The levels agree with each other.** The land is uniform, so level 1 should repeat level 0.
  Over the run the two differ by at most 3.2 W/m² in H, 2.5 W/m² in LE and 0.3 K in the skin
  with Noah-MP; with force-restore by 6.4 W/m² in H, 0.1 K in the skin, and up to 25 W/m² in
  LE, which fluctuates on both levels after the second hour. Neither drifts apart.

Level 1 runs Noah-MP on its own land file (the run prints "Noah-MP at level 1: runs the land
model on its own setup file") and its radiation from its own sweep.
