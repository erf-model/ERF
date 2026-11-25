# ERF domain advisor for HRRR-driven runs

```bash
python3 Exec/Tools/erfdomainadvisor.py --lat 34.07 --lon -118.55 --phenomenon wildfire \
    --wind-from 45 --start 2025-01-07_12 --outdir palisades --plot
```

## Purpose

A run started from HRRR through WPS metgrid needs a domain, a base cell size and a stack of
refined levels.  These are often picked around the point of interest without checking that
the domain is large enough or that the terrain is gentle enough for a terrain-following grid.
This tool takes a point and the phenomenon to be studied, and suggests:

* the outer (level-0) domain and its base cell size;
* the number of refined levels, their cell sizes and their boxes, as WPS nests;
* the vertical grid, the time step and a rough cost;
* with `--fire` (on by default for `wildfire`), the fire grid on the finest level.

It checks:

* the terrain slope on every level;
* the terrain mismatch along the lateral relaxation zone;
* the distance to the HRRR grid edge;
* which levels fall in the turbulence gray zone.

It writes `report.txt`, a `namelist.wps` (d01 plus one nest per refined level), a starting ERF
`inputs` file and `plan.json`.  With `--plot` it also writes `domain.png` and `slopes.png`.

Every number in the tool is a best guess from common downscaling practice, not a value taken
from a paper.  The report says so, and the inputs file is a starting point the user edits.

It sits beside the other two input tools:

* `erfcreateinput.py` writes a deck for a configuration the user has already chosen.
* `check_erf_inputs.py` checks any deck before a run.  Run it on the inputs file this tool
  writes.

## Files

| File | Content |
|---|---|
| `../erfdomainadvisor.py` | Command line; writes the files |
| `phenomena.py` | The phenomena and their best guesses |
| `plan.py` | Level chain, boxes, nests, fire grid, vertical grid, time step, cost; range checks |
| `projection.py` | Lambert conformal projection, the HRRR CONUS grid, the ERF frame |
| `terrain.py` | Elevation sources (AWS Terrain Tiles, ESRI ASCII grid), PNG decoder, slope checks |
| `writers.py` | Report, `namelist.wps`, inputs file, `plan.json` |
| `plot.py` | `--plot`: `domain.png` and `slopes.png` (needs matplotlib) |
| `selftest.py` | `--self-test`, offline |

Python 3.8 or later and numpy are needed.  matplotlib is needed for `--plot` only.

## How to run

```bash
python3 Exec/Tools/erfdomainadvisor.py --list
```

```bash
python3 Exec/Tools/erfdomainadvisor.py --lat 40.0 --lon -105.3 --phenomenon downslope_windstorm --fire --dem none
```

```bash
python3 Exec/Tools/erfdomainadvisor.py --self-test
```

The workflow after the tool:

1. Run WPS (geogrid, ungrib with `Vtable.raphrrr`, metgrid) with the written `namelist.wps` on
   HRRR pressure-level files.  Set `geog_data_path`.
2. Check the deck: `python3 Exec/Tools/check_erf_inputs.py <outdir>/inputs`.
3. Run ERF with the inputs file in the directory of the met_em files.

Main options.  All have range checks that stop the tool with a message naming the option; no
value is clamped.

| Option | Default | Effect |
|---|---|---|
| `--lat`, `--lon` | required | point of interest [deg]; latitude -80 to 80, longitude -180 to 360 |
| `--phenomenon` | required | one of the table below |
| `--target-dx` | phenomenon | finest cell size [m], 1 to 3000 |
| `--ratio` | 3 | refinement ratio, 2 or 3 |
| `--max-base-dx` | 1000 | largest base cell [m], 300 to 3000 |
| `--outer LX LY` | phenomenon | smallest outer domain [km], 1 to 5000 |
| `--finest LX LY` | phenomenon | finest box [km], 1 to 5000 |
| `--fetch` | phenomenon | upwind fetch added to the finest box [km], 0 to 200 |
| `--wind-from` | none | direction the wind comes from [deg], 0 to below 360; puts the fetch upwind |
| `--top`, `--dz0`, `--hours` | phenomenon | model top [km] 1 to 40, first cell [m] 0.5 to 200, run length [h] 0.25 to 240 |
| `--start` | none | first time, a whole hour, `YYYY-MM-DD_HH` (or `YYYY-MM-DD_HH:00:00`); names the met_em files |
| `--fire` / `--no-fire` | phenomenon | write an `erf.fire.*` block, fire grid on the finest level; give one or neither |
| `--fire-dx` | 10 | largest fire cell [m], 0.5 to 100 |
| `--real-width` | 5 | `erf.real_width` [cells], 1 to 50 |
| `--nest-margin` | 12 | smallest distance of a box to its parent's edge [parent cells], 4 to 100 |
| `--hrrr-buffer` | 5 | HRRR cells between the relaxation zone and the first box, and to the HRRR edge, 0 to 100 |
| `--max-stretch` | 1.08 | largest vertical stretching ratio, above 1 to 1.3 |
| `--cost-per-cell-step` | 2e-6 | core-seconds per cell and step for the cost estimate, 1e-9 to 1e-3 |
| `--dem` | `terrarium` | `terrarium`, an ESRI ASCII grid in degrees, or `none` |
| `--slope-caution`, `--slope-limit` | 0.3, 0.6 | slopes that get a note and a warning |
| `--mismatch` | 150 | boundary terrain mismatch that gets a warning [m] |
| `--outdir`, `--force`, `--plot` | | where to write; overwrite; write the images |

## Rules (best guesses)

* **Base cell.**  The finest cell times `ratio**L`, with `L` the largest number of levels that
  keeps the base at or below `--max-base-dx` (1000 m, a third of the HRRR spacing).  The
  targets are 900 m divided by powers of 3, so the default chain starts at 900 m.  More than 6
  refined levels are refused.
* **Boxes.**  The finest box is the phenomenon's box around the point plus the fetch: all of
  it upwind with `--wind-from`, half on each side without.  Each coarser box adds
  `--nest-margin` cells of its own size.
* **Outer domain.**  At least the phenomenon's outer size, shifted upwind by 15 % of its
  width.  It keeps the first refined box `erf.real_width` base cells plus `--hrrr-buffer` HRRR
  cells plus `--nest-margin` base cells inside its edges.
* **Nests.**  Boxes are snapped outward onto their parent's cells, so each is a WPS nest with
  `e_we - 1` divisible by the ratio.  ERF reads the place of each refined level from its
  met_em file, so the inputs file needs no refinement boxes.
* **Fire grid.**  `erf.fire.grid_ratio` is the smallest integer that makes the fire cell at most
  `--fire-dx`.  ERF needs every x,y box length on the fire level to be a multiple of it.  The
  finest box is grown to a multiple of a box length near 64 cells (a multiple of the grid ratio
  and the refinement ratio).  The deck sets `amr.max_grid_size_x/y` to that length there, so
  AMReX cuts the level into equal boxes.  The parents are padded so that the nest margin holds.
* **Vertical grid.**  The smallest `nz`, a multiple of 8, with a stretching ratio at or below
  `--max-stretch` and a top cell at or below the phenomenon's largest cell.  The deck gives
  `erf.initial_dz` only; ERF solves the ratio.  When the first cell is already as large as
  the top cell allows, the grid is uniform, the deck omits `erf.initial_dz` and the report
  says so.
* **Turbulence.**  `dx >= zi`: PBL scheme (MYNN 2.5).  `dx < 0.2 zi`: LES (Smagorinsky, no PBL
  scheme).  In between (`0.2 zi <= dx < zi`): gray zone, reported.  Set per level in `erf.pbl_type` and
  `erf.les_type`.
* **Time step.**  The largest step that divides 900 s and stays at or below both 6 s per km
  (WRF practice) and the limit ERF applies to a fixed step (`fixed_dt_cfl_limit`, on by
  default): 0.8 dx / (m (sqrt(2) u_max + 2 c / 3)), with c = 370 m/s (air at about 340 K, a margin over
  near-surface air) and m the largest map factor in the domain.  Over a fire the air is
  hotter and ERF may lower the step there.  Fine levels subcycle by the ratio.  `erf.fixed_dt` keeps the
  15-minute output on a step; the deck says how to switch to `erf.cfl`.
* **Terrain.**  Cell means on each level's grid, with water clipped at 0 m.  The slope of a cell
  is taken from one-cell drops (the larger of the two neighbours in x and in y), as the
  terrain-following grid sees it.  The slope note is at 0.3, the warning at 0.6.  A level
  with more than 4 million cells is sampled every few cells, and the report says so.  The
  report also gives ERF's terrain slope factor alpha = slope dx / dz1, with dz1 = dz0 (top - h) / top
  the first cell over the highest ground h.  ERF reports alpha in its log when
  `erf.diffusive_dt_check = true`.  Above 1 the explicit terrain-metric diffusion limits the
  step, and `erf.implicit_terrain_metric` is suggested.  The report warns when the highest
  ground leaves less than the boundary layer plus the damping layer below the top.  When the finest cell is below 900 m, the
  report notes that WPS should use terrain finer than its default GMTED2010 30" data
  (`geog_data_res`).
* **Terrain smoothing.**  The deck sets `erf.terrain_smoothing = 0`.  On the Metgrid path ERF
  stops on a refined level when it is not 0, because each level's mesh comes from its own
  met_em file.  Isolated spikes in the tiles (more than 200 m from the 3 x 3
  median) are replaced by the median.  The boundary mismatch compares level-0 terrain with its
  3 km average over the relaxation zone.
* **Cost.**  Cells times steps times 2e-6 core-seconds, and 2000 bytes per cell.  Both are guesses.

## Phenomena

The values are best guesses.  `--list` prints the same table; every value can be overridden.

| Phenomenon | Finest dx [m] | Outer [km] | Finest box [km] | Fetch [km] | Top [km] | dz0 [m] | Hours | zi [m] | Moisture | Fire |
|---|---|---|---|---|---|---|---|---|---|---|
| `convective_bl_les` | 100 | 60 x 60 | 30 x 30 | 15 | 6 | 20 | 10 | 1500 | Kessler | - |
| `deep_convection` | 300 | 400 x 400 | 150 x 150 | 0 | 20 | 40 | 12 | 1500 | Morrison | - |
| `downburst_gust` | 33.33 | 150 x 150 | 15 x 15 | 0 | 15 | 5 | 6 | 1500 | Morrison | - |
| `downslope_windstorm` | 100 | 300 x 200 | 40 x 40 | 10 | 20 | 20 | 12 | 1000 | Kessler | - |
| `fog_stratus` | 300 | 150 x 150 | 60 x 60 | 0 | 4 | 5 | 18 | 300 | Morrison | - |
| `haboob_dust` | 300 | 300 x 300 | 120 x 120 | 0 | 12 | 20 | 12 | 2000 | Morrison | - |
| `hill_flow` | 33.33 | 80 x 80 | 8 x 8 | 5 | 6 | 3 | 6 | 1000 | Kessler | - |
| `hurricane_landfall` | 1000 | 1000 x 1000 | 400 x 400 | 0 | 20 | 50 | 24 | 1000 | Morrison | - |
| `lake_effect_snow` | 300 | 400 x 300 | 150 x 100 | 0 | 12 | 30 | 24 | 1000 | Morrison | - |
| `mountain_waves` | 300 | 400 x 300 | 150 x 100 | 0 | 25 | 50 | 12 | 1000 | Kessler | - |
| `sea_breeze` | 300 | 200 x 200 | 80 x 80 | 0 | 6 | 20 | 12 | 1000 | Kessler | - |
| `urban_canopy` | 3.704 | 100 x 100 | 3 x 3 | 3 | 4 | 2 | 3 | 1000 | Kessler | - |
| `valley_cold_pool` | 33.33 | 80 x 80 | 15 x 15 | 0 | 8 | 5 | 18 | 300 | Kessler | - |
| `wildfire` | 33.33 | 150 x 150 | 15 x 15 | 10 | 10 | 10 | 12 | 1500 | Kessler | yes |
| `wind_farm` | 33.33 | 100 x 100 | 25 x 20 | 15 | 6 | 10 | 12 | 1000 | Kessler | - |

## Expected results

For the wildfire example at the top (34.07 N, 118.55 W, wind from 45 deg):

| Level | dx [m] | Cells | Regime |
|---|---|---|---|
| 0 | 900 | 168 x 168 x 80 | gray zone (PBL scheme) |
| 1 | 300 | 120 x 108 x 80 | gray zone |
| 2 | 100 | 279 x 246 x 80 | LES |
| 3 | 33.3 | 720 x 660 x 80 | LES, fire grid ratio 4 (8.3 m) |

The base time step is 2 s.  The steepest one-cell slope is 0.71 at 900 m and 1.89 at 33 m,
so the report warns about every level; the slope factor alpha is 90 on level 0 and 6.8 on
level 3, so it suggests `erf.implicit_terrain_metric`.  The boundary mismatch is 227 m.  The checker reports only the met_em files, which do not exist
yet, and advice to prefer `erf.cfl`.

## CTest

`DomainAdvisor_SelfTest` (label `unit`) runs `--self-test`.  It is offline and needs no ERF
build.  Without numpy, the 23 tests that do not need it run and the 9 terrain tests are
skipped.  The tests compare against closed forms or published values:

| Test | Reference |
|---|---|
| HRRR corners | first point 21.138123 N, 122.719528 W; far corner 47.842195 N, 60.917193 W |
| Lambert round trip | inverse(forward(p)) = p; map factor 1 on the standard parallel |
| Level chain | 33.3 m by 3 from 900 m (3 levels); 100 m by 2 from 800 m |
| Vertical grid | dz0 (q^nz - 1) / (q - 1) = top; nz a multiple of 8 and the smallest |
| Nests and ERF box | the box ERF builds from the met_em attributes holds the point; sizes and margins |
| `namelist.wps` | `e_we`, `e_sn`, `i_parent_start` match the plan |
| Fire grid | grid ratio 4 for 33.3 m and 10 m; fire level lengths multiples of the box length |
| Gaussian hill | steepest slope H / s exp(-1/2) from an ESRI ASCII grid |
| Slope warnings | a ramp of slope 0.8 is warned about on every level; flat ground is not |
| Despike, mismatch | one spike replaced; a 500 m cell gives 500 - 500 / 9 m |
| PNG decoder | round trip through filters 0 to 4, and a Paeth tie |
| Deck | the written deck has no checker error or warning apart from the met_em files |
| Deck cell size | `geometry.prob_hi / amr.n_cell` equals the met_em DX to within ERF's 1e-3 m |
| Deck smoothing | `erf.terrain_smoothing = 0` with refined levels |
| Time step | 900 m and 40 m/s give 2 s, below ERF's fixed-step limit (2.37 s) |
| Vertical edge cases | a uniform grid, and a column of more than 1024 cells without overflow |
| Fire grid, ratio 2 | the padded parents keep the nest margin |
| One-cell slope | a 100 m step over 10 m cells gives 10 (central differences give 5) |
| Slope factor | alpha = slope dx / dz1 with the first cell squeezed over 100 m of ground; above 1 is warned about; a 5 km plateau under a 6 km top is warned about |
| Damaged tiles | a cut or corrupt PNG gives a terrain error, not a traceback |
| Regime | dx = 0.2 zi is gray zone, just below is LES |
| Gray-zone text | the warnings state `0.2 zi <= dx < zi`, as the code decides |
| Input ranges | each option outside its range aborts with its name |
| ERF stretching solve | 48 cells of 0.5 m to 4 km need more than ERF's 50 Newton steps; the plan adds cells until ERF's iteration converges; a uniform grid meets the top-cell limit |
| Outside HRRR | a domain past the HRRR edge is refused |
| Antimeridian tiles | a point at 180 deg fetches at most 4 tiles |
| Sampled mismatch | the boundary mismatch uses the sampling spacing on a sampled level |
| Box length | the box length on a refined level is a multiple of the ratio |
| `--outdir` a file | exit code 2 and a message, no traceback |
| Fire flags | `--fire` with `--no-fire` is refused |
| Start time | the hour is needed; minutes or seconds other than 0 are refused |
| Shallow top | a top below the boundary layer plus the damping layer is warned about |

Each test was seen to fail against a one-line mutant of the code it covers: 50 mutants, all
caught.  One more mutant (Paeth `pa < pb` for `pa <= pb`) gives the same output, because
`pa = pb` with a different left and up value forces `pc = 0`.

## Data and acknowledgements

**Elevation: AWS Terrain Tiles (terrarium encoding).**  Made by Mapzen and the Tilezen project
(<https://github.com/tilezen/joerd>), hosted by Amazon Web Services as an open data set
(<https://registry.opendata.aws/terrain-tiles/>).  We thank Mapzen, Tilezen and AWS, and the
data providers below.  The tool downloads the tiles it needs into
`~/.cache/erf_domain_advisor/terrarium` at run time.  No tile is stored in this repository.
The attribution the providers require (from
<https://github.com/tilezen/joerd/blob/master/docs/attribution.md>):

* ArcticDEM terrain data DEM(s) were created from DigitalGlobe, Inc., imagery and funded under
  National Science Foundation awards 1043681, 1559691, and 1542736;
* Australia terrain data © Commonwealth of Australia (Geoscience Australia) 2017;
* Austria terrain data © offene Daten Österreichs – Digitales Geländemodell (DGM) Österreich;
* Canada terrain data contains information licensed under the Open Government Licence – Canada;
* Europe terrain data produced using Copernicus data and information funded by the European
  Union - EU-DEM layers;
* Global ETOPO1 terrain data U.S. National Oceanic and Atmospheric Administration;
* Mexico terrain data source: INEGI, Continental relief, 2016;
* New Zealand terrain data Copyright 2011 Crown copyright (c) Land Information New Zealand and
  the New Zealand Government (All rights reserved);
* Norway terrain data © Kartverket;
* United Kingdom terrain data © Environment Agency copyright and/or database right 2015. All
  rights reserved;
* United States 3DEP (formerly NED) and global GMTED2010 and SRTM terrain data courtesy of the
  U.S. Geological Survey.

Images made with `--plot` from these tiles carry the same attribution.

**HRRR grid.**  The tool uses the HRRR CONUS grid layout (Lambert conformal, 3 km,
1799 x 1059 points, standard parallel 38.5 N, meridian 97.5 W).  It does not read HRRR data.
HRRR is run by NOAA (Dowell et al. 2022, *Weather and Forecasting* 37, 1371-1395,
<https://doi.org/10.1175/WAF-D-21-0151.1>; Benjamin et al. 2016, *Monthly Weather Review* 144,
1669-1694, <https://doi.org/10.1175/MWR-D-15-0242.1>).  We thank NOAA for the HRRR system and
its open archive.

## Not included

* Reading HRRR files, or running WPS.
* Fire-engine keys (ignition, fuel map, fuel moisture), wind-farm keys and building keys: they
  belong to their modules and are left to the user.
* Projections other than Lambert conformal.  Near the equator the report suggests Mercator,
  which the user sets by hand.
* Domains for HRRR Alaska, or for GFS outside the HRRR grid: the report warns when the point is
  outside HRRR CONUS.
