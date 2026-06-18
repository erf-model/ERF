# ERF Sep-2024 cyclone validation against NEXRAD / MRMS / ASOS

Compare ERF model output to observations using **event aggregates** (max / total
over the event) so the comparison is robust to the very different cadences
(ERF 30 min, MRMS 2 min, NEXRAD ~5 min, ASOS ~5 min).

| Quantity            | ERF field      | Obs source                        |
|---------------------|----------------|-----------------------------------|
| peak 10 m wind      | `wspd10max`    | ASOS peak gust / 2-min wind       |
| storm total precip  | `rain_accum`   | MRMS `RadarOnly_QPE_01H` summed   |
| event-max reflect.  | `cmpref_max`   | MRMS composite max / NEXRAD comp. |

All ERF fields are **running aggregates updated every step** (cadence-independent;
`wspd10max`/`cmpref_max` from the SurfaceLayer diagnostics, `rain_accum` from the
microphysics). The MOST log-law 10 m wind / 2 m T are in `u10,v10,wspd10,t2,q2`.

## 1. Download observations (no auth)
```
python get_asos.py   <outdir>/asos        # IEM ASOS CSV (14 coastal stations)
python get_mrms.py   <outdir>/mrms --all  # NOAA AWS noaa-mrms-pds (composite refl + QPE)
python get_nexrad.py <outdir>/nexrad --all# Unidata mirror unidata-nexrad-level2 (L2 volumes)
```
NOTE: the NOAA `noaa-nexrad-level2` bucket blocks anonymous listing; we use the
`unidata-nexrad-level2` mirror (same layout, anonymous list+GET).

## 2. Event aggregates
```
python aggregate_event.py asos <asos_dir>          -> asos_event.csv (per-station peaks)
```

## 3. Put everything on the ERF grid and compare
Run from the ERF run directory (needs plt2d_* with lat_m/lon_m + the diagnostics):
```
python compare_obs.py erf-grid    "plt2d_sep[0-9]*" "plt_sep[0-9]*"     -> erf_grid.npz
python compare_obs.py mrms-grid   <mrms_dir>   erf_grid.npz             -> mrms_on_erf.npz
python compare_obs.py nexrad-grid <nexrad_dir> erf_grid.npz            -> nexrad_on_erf.npz
python compare_obs.py asos-pts    asos_event.csv erf_grid.npz          -> asos_compare.csv
python compare_obs.py plot erf_grid.npz mrms_on_erf.npz nexrad_on_erf.npz asos_compare.csv
        -> compare_reflmax.png, compare_asos_wind.png
```
Regridding obs -> ERF cells is nearest-neighbour (scipy.cKDTree); a mild
coarsening from MRMS 1 km / NEXRAD polar to ERF 3/9 km, appropriate for maxima.
ERF boundary relaxation ring (9 cells) is masked via the `interior` array.

## Caveats
- ASOS *gust* is a higher statistic than a sustained 10 m max; `wspd10max` is the
  sustained 10 m peak. Compare to ASOS 2-min wind for like-for-like, gust as upper ref.
- A freshly cold-started run has low `wspd10max` until the storm spins up; only the
  full-event aggregate is meaningful.
- NEXRAD composite here is column-max over sweeps gridded to ERF cells (approx,
  not a full 3D objective analysis) — fine for the event-max swath comparison.
