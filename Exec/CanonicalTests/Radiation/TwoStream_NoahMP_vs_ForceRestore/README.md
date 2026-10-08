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
| Canopy | Jarvis stomata (`CANOPY_STOMATAL_RESISTANCE_OPTION = 2`) | Noah's big-leaf Jarvis resistance on Noah-MP's option-2 parameters (RS, RGL, HS, TOPT, RSMAX), LAI 2.23 (the table's at 15:00 UTC on 5 August, day 217.625 counted from 0 as Noah-MP does); Noah-MP itself applies them per sunlit and shaded leaf |
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
vegetation) on the same soil at its wilting point (0.12 m3/m3) under both models. The pore
air there is nearly dry (relative humidity about 0.005), so neither evaporates and what
remains is the dry energy balance. The force-restore case has no vegetation type: its bare
soil evaporates from its pore air through Noah-MP's soil resistance, its roughness is Z0SOIL
(0.002 m), and it takes the albedo and emissivity Noah-MP gives this soil (0.197 and 0.97).
Its deck is `inputs_force_restore_barren`.

The **moist barren** variant (`--barren-moist`) is the same bare land at 0.25 m3/m3, where
the pore air is nearly saturated and the soil resistance sets the evaporation of both models. Noah-MP's bare-soil albedo falls
with the top layer's water content (`GroundAlbedoMod`): 0.145 rising to 0.151 through the day
here as the top soil dries, so the deck `inputs_force_restore_barren_moist` takes the
shortwave-weighted 0.148. The force-restore decks read the settings
they share from `inputs_force_restore_common`.

## Running it

Needs an executable built with `-DERF_ENABLE_NOAHMP=ON` (and so a parallel NetCDF), `ncgen`
on the `PATH`, and Python 3 with numpy, matplotlib and yt:

```
./run_comparison.sh /path/to/erf_exec [nranks] [--barren | --barren-moist]
```

It writes the land files, runs both cases under `runs/` (`runs_barren/`, `runs_barren_moist/`), and calls
`compare.py`, which writes `comparison.csv` (level means every 10 minutes) and
`comparison.png` there.

## Results

From the release build on 2 MPI ranks (about 12 minutes per case). Level-0 means, Noah-MP
vs force-restore:

| Grassland | 18:50 UTC (3.8 h, Noah-MP's peak H) | 21:00 UTC (6 h) |
|---|---|---|
| Absorbed shortwave [W/m²] | 780.6 vs 771.6 | 633.7 vs 636.8 |
| Net longwave, down [W/m²] | -92.7 vs -82.2 | -89.0 vs -80.4 |
| Skin temperature [K] | 307.8 vs 306.1 | 308.0 vs 306.6 |
| H [W/m²] | 115.9 vs 128.7 | 79.4 vs 77.5 |
| LE [W/m²] | 377.9 vs 397.7 | 348.9 vs 354.4 |
| PBL height [m] | 958 vs 944 | 1152 vs 1156 |

| Barren | 19:10 UTC (4.2 h, Noah-MP's peak H) | 21:00 UTC (6 h) |
|---|---|---|
| Absorbed shortwave [W/m²] | 779.0 vs 778.4 | 647.4 vs 647.3 |
| Net longwave, down [W/m²] | -208.7 vs -197.0 | -196.5 vs -186.4 |
| Skin temperature [K] | 326.6 vs 324.7 | 326.4 vs 324.7 |
| H [W/m²] | 364.0 vs 345.1 | 308.8 vs 293.0 |
| LE [W/m²] | -4.1 vs -3.4 | -3.8 vs -3.2 |
| PBL height [m] | 1477 vs 1406 | 1809 vs 1713 |

At the wilting point both bare soils take up a little vapour (LE of -3 to -4 W/m²): the
pore air is nearly dry (relative humidity about 0.005), so the bare-soil evaporation follows
Noah-MP's pore-air humidity, not only its resistance.

| Moist barren | 19:10 UTC (4.2 h, Noah-MP's peak H) | 21:00 UTC (6 h) |
|---|---|---|
| Absorbed shortwave [W/m²] | 825.1 vs 825.9 | 684.4 vs 686.8 |
| Net longwave, down [W/m²] | -157.8 vs -146.4 | -150.1 vs -138.9 |
| Skin temperature [K] | 318.7 vs 316.9 | 318.7 vs 316.7 |
| H [W/m²] | 233.5 vs 212.3 | 195.4 vs 174.0 |
| LE [W/m²] | 205.0 vs 205.9 | 182.9 vs 184.8 |
| PBL height [m] | 1201 vs 1165 | 1494 vs 1410 |

At 0.25 m3/m3 the pore air is nearly saturated and the bare-soil resistances set the
evaporation: the two LE agree within 1 %. The force-restore skin is 2 K cooler and its H
about 10 % lower, as on the other surfaces.

On all three surfaces the force-restore balance comes within about 2 K of Noah-MP's skin,
10-15 % of its H and LE, and about 5 % of its boundary-layer depth. Two choices decide that, and both now
come from Noah-MP's tables:

- **The evaporation.** The tables above compare against Noah-MP with Jarvis stomata
  (`CANOPY_STOMATAL_RESISTANCE_OPTION = 2`), the form the balance has. Noah-MP's default,
  Ball-Berry (option 1), transpires much less on this day. The grass LE at 3.8 h, against
  each reference:

  | Grass LE at 3.8 h [W/m²] | |
  |---|---|
  | Force-restore, no moisture coupling (fixed surface mixing ratio)* | 322 |
  | Force-restore, soil-water factor alone* | 474 |
  | Force-restore, resistances (this case) | 398 |
  | Noah-MP, Jarvis (option 2, this case) | 378 |
  | Noah-MP, Ball-Berry (option 1, the default) | 166.5 |

  \* From runs made while the coupling was developed, with the 0.1 m roughness of the earlier
  decks.

  Without the coupling, ERF requires `erf.most.surf_moist` to equal the sounding's surface
  value, so the surface evaporates from a fixed mixing ratio. The soil-water factor alone
  scales the evaporation from q_sat at the skin, which is well above that fixed value. The
  canopy and bare-soil resistances bring the balance within about 5 % of Noah-MP with the same
  stomata. Against default Noah-MP (Ball-Berry: LE 166.5, H 282.5 W/m² and a skin of
  313.4 K) its LE is still 2.4 times higher. The agreement is with Noah-MP using the same
  canopy model, not with every Noah-MP configuration.
- **The roughness.** Over bare soil, with the 0.1 m of the earlier decks, the balance's
  surface sheds its heat far more easily than Noah-MP's 0.002 m soil: a skin 10.7 K cooler and
  30 % more H at the peak. With the tables' roughness the difference is about 2 K and 5 %.

The balance has the Jarvis form only, which is why this case runs Noah-MP with option 2.

Both models now count the day of the year from 0, so they read the same table LAI, 2.23 at the
start. The driver used to count from 1 and read 2.16 here; the flux and skin differences quoted
below were measured then, and are expected to shift a little the next time the case is run.

Within each model, level 1 tracks level 0 (the land is uniform). On grass the two differ by
at most about 5 W/m² in H and LE and 0.2 K in the skin. On bare soil the skins stay within 0.5 K,
but H differs by up to 9 W/m² (Noah-MP) and 36 W/m² (force-restore) at moments in the
afternoon; by the end of the run every difference is under 5 W/m², so the levels do not
drift apart. On moist bare soil Noah-MP's levels differ by up to 14 W/m² in H (8 W/m² at the
end) and the force-restore levels by under 4 W/m².
