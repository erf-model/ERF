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

| | `inputs_noahmp` | `inputs_force_restore` (on `inputs_force_restore_common`) |
|---|---|---|
| Surface model | Noah-MP on both levels | force-restore balance on both levels |
| Land | grassland (MODIS 10) on 80 % of the surface (`SHDMAX`), silty clay loam (STAS 8), 290 K soil at 0.25 m3/m3; level 1 has its own land setup file | the same, from Noah-MP's tables: `seb_soil_type = 8`, `seb_vegetation_type = 10`, `seb_vegetation_fraction = 0.8`, soil water 0.25 in the top 0.1 m |
| Canopy | Jarvis stomata (`CANOPY_STOMATAL_RESISTANCE_OPTION = 2`) | Noah's big-leaf Jarvis resistance on Noah-MP's option-2 parameters (RS, RGL, HS, TOPT, RSMAX), LAI 2.2 (the table's for 5 August); Noah-MP itself applies them per sunlit and shaded leaf |
| Bare soil | Noah-MP's soil resistance (option 1) | the same resistance |
| Evaporation | Noah-MP's | the surface layer's, from beta q_sat(T_s) + (1 - beta) q_air (`seb_surface_layer_uses_moisture`); the balance drains its soil water by it |
| Roughness | Noah-MP's own (it supplies u*, theta*, q*) | 0.0964 m, 0.8 Z0MVT + 0.2 Z0SOIL from the tables |
| Surface albedo and emissivity | Noah-MP's: emissivity 0.994; albedo 0.229 at 15:10 UTC falling to 0.201 at 18:00 | 0.21 (about Noah-MP's daytime mean) and 0.994 |
| Ground heat | Noah-MP's soil column | restore to 290 K (the soil's) over a day, with the surface heat capacity of the soil (1.16e5 J/m2/K, from its heat capacity and conductivity) |
| Land step | 60 s (`NOAH_TIMESTEP`) | every ERF step |
| H and LE | Noah-MP's, applied through the surface layer | the surface layer's, which the balance removes; the surface layer takes its temperature and moisture from the balance (`seb_surface_layer_uses_skin`, `seb_surface_layer_uses_moisture`) |

Noah-MP runs on both levels because `namelist.erf` names a land setup file for each
(`ERF_SETUP_FILE_01` and `_02`). `make_land_files.py` writes them; level 1's places itself
in level 0 with the WRF attributes `I_PARENT_START = J_PARENT_START = 3` and
`PARENT_GRID_RATIO = 2`.

The **barren** variant (`run_comparison.sh ... --barren`) puts bare land (MODIS 16, no
vegetation) on the same soil at its wilting point (0.12 m3/m3) under both models, so that
neither evaporates and what remains is the dry energy balance. The force-restore case then
has no vegetation type (its roughness is Z0SOIL, 0.002 m), and takes the albedo and
emissivity Noah-MP gives this soil (0.197 and 0.97); its deck is `inputs_force_restore_barren`.
Both force-restore decks read the settings they share from `inputs_force_restore_common`.

## Running it

Needs an executable built with `-DERF_ENABLE_NOAHMP=ON` (and so a parallel NetCDF), `ncgen`
on the `PATH`, and Python 3 with numpy, matplotlib and yt:

```
./run_comparison.sh /path/to/erf_exec [nranks] [--barren]
```

It writes the land files, runs both cases under `runs/` (`runs_barren/`), and calls
`compare.py`, which writes `comparison.csv` (level means every 10 minutes) and
`comparison.png` there.

## Results

From the release build on 2 MPI ranks (about 12 minutes per case). Level-0 means, Noah-MP
vs force-restore:

| Grassland | 18:50 UTC (3.8 h, Noah-MP's peak H) | 21:00 UTC (6 h) |
|---|---|---|
| Absorbed shortwave [W/m²] | 780.6 vs 771.5 | 633.7 vs 636.8 |
| Net longwave, down [W/m²] | -92.7 vs -82.5 | -89.0 vs -80.5 |
| Skin temperature [K] | 307.8 vs 306.1 | 308.0 vs 306.6 |
| H [W/m²] | 115.9 vs 130.1 | 79.4 vs 78.0 |
| LE [W/m²] | 377.9 vs 395.7 | 348.9 vs 353.9 |
| PBL height [m] | 958 vs 945 | 1152 vs 1159 |

| Barren | 19:10 UTC (4.2 h, Noah-MP's peak H) | 21:00 UTC (6 h) |
|---|---|---|
| Absorbed shortwave [W/m²] | 779.0 vs 778.4 | 647.4 vs 647.3 |
| Net longwave, down [W/m²] | -208.7 vs -196.3 | -196.5 vs -185.9 |
| Skin temperature [K] | 326.6 vs 324.6 | 326.4 vs 324.6 |
| H [W/m²] | 364.0 vs 343.2 | 308.8 vs 288.4 |
| LE [W/m²] | -4.1 vs 0.0 | -3.8 vs 0.0 |
| PBL height [m] | 1477 vs 1401 | 1809 vs 1720 |

On both surfaces the force-restore balance comes within about 2 K of Noah-MP's skin, 10-15 %
of its H and LE, and 5 % of its boundary-layer depth. Two choices decide that, and both now
come from Noah-MP's tables:

- **The evaporation.** Without the moisture coupling the balance's surface evaporates from a
  fixed mixing ratio (ERF requires `erf.most.surf_moist` to equal the sounding's surface
  value) and so behaves as a wet surface: at 3.8 h, LE 322 against Noah-MP's 166 W/m²
  (Ball-Berry) and a skin 9 K cooler. The soil-water factor alone makes it worse (LE 474): it
  scales evaporation from q_sat at the skin, which is well above the fixed value. The canopy
  and bare-soil resistances bring it to Noah-MP's.
- **The roughness.** Over bare soil, with the 0.1 m of the earlier decks, the balance's
  surface sheds its heat far more easily than Noah-MP's 0.002 m soil: a skin 10.7 K cooler and
  30 % more H at the peak. With the tables' roughness the difference is 2 K and 6 %.

Noah-MP's default stomata, Ball-Berry (`CANOPY_STOMATAL_RESISTANCE_OPTION = 1`), transpire
much less here: at 3.8 h, LE 166.5 and H 282.5 W/m² with a skin of 313.4 K. The balance has
the Jarvis form only, which is why this case runs Noah-MP with option 2.

Within each model, level 1 tracks level 0 (the land is uniform). On grass the two differ by
at most 5 W/m² in H and LE and 0.2 K in the skin. On bare soil the skins stay within 0.5 K,
but H differs by up to 9 W/m² (Noah-MP) and 44 W/m² (force-restore) at moments in the
afternoon; by the end of the run every difference is under 5 W/m², so the levels do not
drift apart.
