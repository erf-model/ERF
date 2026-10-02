# Noah-MP from an idealized case

This case runs the Noah-MP land-surface model under an idealized ERF atmosphere: an
`input_sounding` initialization on a 4 x 4 x 32 periodic column, with no radiation model.
It is the only case that builds and runs Noah-MP without a WRF or WPS initialization.

The atmosphere does not have to come from WRF, but Noah-MP's driver (WRF's, with ERF
patches) always reads its land state from a wrfinput-format NetCDF file. Here that file
describes one grassland patch (vegetation type 10) on silty clay loam (soil type 8), with
four soil layers.

## Build

Noah-MP needs NetCDF built with parallel I/O (see
[Coupling to Noah-MP](https://erf.readthedocs.io/en/latest/CouplingToNoahMP.html)):

```
cmake -DERF_ENABLE_MPI=ON -DERF_ENABLE_NETCDF=ON -DERF_ENABLE_NOAHMP=ON ...
```

or, with GNU Make from `Exec/` (NetCDF, including netcdf-fortran, found through
`pkg-config`):

```
make USE_MPI=TRUE USE_NETCDF=TRUE USE_NOAHMP=TRUE
```

## Run

Everything the run needs is in this directory, so it runs in place:

```
mpiexec -n 1 /path/to/erf_exec inputs_noahmp_ideal
```

Noah-MP's driver reads three files from the run directory:

- `namelist.erf`: the driver options. It names `wrfinput_d01` as the land setup file
  (`ERF_SETUP_FILE_01`), and its `soil_thick_input` matches `DZS` in that file. It must not
  set `ZLVL`, since ERF passes the reference height for every column.
- `wrfinput_d01`: the land setup file, a 4 KB classic NetCDF file. It is generated from
  `wrfinput_ideal.cdl`, the readable, commented source of the same data. After changing the
  CDL, regenerate it with

  ```
  ncgen -o wrfinput_d01 wrfinput_ideal.cdl
  ```

- `NoahmpTable.TBL`: Noah-MP's parameter table, a copy of
  `Submodules/Noah-MP/parameters/NoahmpTable.TBL` at the submodule commit this case was
  added with. Refresh it from the submodule if Noah-MP's table format changes.

## What to expect

The run takes two 1 s ERF steps. Noah-MP fires once, on the first step, and advances one
3600 s land step (`NOAH_TIMESTEP`). `plt00002` holds the atmosphere at the last step;
`plt2d00000` and `plt2d00002` hold Noah-MP's outputs `t_sfc`, `sav`, `sag`,
`sensible_heat_flux`, `grdflx` and `fira`.

No radiation model is set, so nothing supplies the downwelling shortwave, downwelling
longwave or solar zenith angle Noah-MP reads. ERF prints a warning at start-up and another
at the first land step, and passes zero for all three. In the output:

- the absorbed shortwave, `sav` and `sag`, is exactly zero;
- the surface cools from the 300 K in the setup file to about 252 K. With no downwelling
  longwave it radiates to a 0 K sky, and it settles where its own emission is balanced by
  heat from below and above:

  | Term at the last step (positive away from the surface) | W/m² |
  | --- | --- |
  | net longwave, `fira` | +228.4 |
  | ground heat flux, `grdflx` (heat conducted up from the warmer soil) | −151.2 |
  | sensible heat, `sensible_heat_flux` (heat from the warmer air) | −78.2 |
  | latent heat (evaporation), Noah-MP `LH` | +1.0 |
  | absorbed shortwave, `sav` + `sag` | 0 |

  Noah-MP's land output, `lnd00002/Level_0.nc`, closes this budget (`FIRAXY`, `HFX`, `LH`,
  `GRDFLX`) to 0.02 W/m², so evaporation plays almost no part even though the air is dry.
  The latent heat is not in `plt2d00002`: ERF fills its `latent_heat_flux` there only when a
  moisture model is set. This is the correct response to the configuration, not a Noah-MP
  defect.

A run that stops before `plt00002` and `plt2d00002` are written has failed, whatever its
exit status. Noah-MP's own physics checks end with a Fortran `STOP`, which ERF reports as a
failure.

## With two-stream radiation: `inputs_noahmp_twostream`

The same land setup and atmosphere, with `erf.radiation_model = TwoStream` supplying
Noah-MP's radiative forcing. The time is 18:00 UTC on 2024-08-05 at the setup file's site
(`erf.rad_cons_lat = 40`, `erf.rad_cons_lon = -100`), about 11:20 local solar time, so the
sun is well up. The two-stream optical depths are set for this 1 km column (0.1 in total for
the shortwave, 1.6 for the longwave) rather than left at their whole-atmosphere defaults.
It runs in place like the case above:

```
mpiexec -n 1 /path/to/erf_exec inputs_noahmp_twostream
```

Every step the two-stream model writes the downwelling shortwave and longwave at the surface
and the cosine of the solar zenith angle into Noah-MP's inputs. `plt2d00002` carries them
as `sw_flux_dn`, `lw_flux_dn` and `cos_zenith_angle` (a run without the coupling writes
`-999` there), next to Noah-MP's outputs:

| Field in `plt2d00002` | Value |
| --- | --- |
| `cos_zenith_angle` | 0.9077 |
| `sw_flux_dn` (SWDOWN) | 1075.4 W/m² |
| `lw_flux_dn` (GLW) | 354.4 W/m² |
| `t_sfc` | 309.857 K (252.306 K without radiation) |
| `sav` + `sag` | 441.7 + 417.4 W/m² (0 without radiation) |

The surface now warms instead of cooling, and Noah-MP's land output `lnd00002/Level_0.nc`
closes its energy budget (positive away from the surface, except the absorbed shortwave):

| Term | W/m² |
| --- | --- |
| absorbed shortwave, `FSAXY` (= `SAVXY` + `SAGXY`) | 859.18 |
| net longwave, `FIRAXY` | 167.21 |
| sensible heat, `HFX` | 249.38 |
| latent heat, `LH` | 147.97 |
| ground heat flux into the soil, `GRDFLX` | 294.63 |
| residual, `FSAXY` − (`FIRAXY` + `HFX` + `LH` + `GRDFLX`) | −0.012 |

The residual is the same size as in the case without radiation (0.016 W/m²).

In the other direction the two-stream model takes Noah-MP's broadband `albedo` (0.201 here)
as its surface albedo, so it reflects the shortwave Noah-MP reflects: with
`erf.radiation.diag_csv_enable = true` the radiation CSV reports `SW_surface` = 859.19 W/m²
absorbed at the ground at step 1, Noah-MP's `sav` + `sag` = 859.18 W/m². (The visible
direct-beam albedo `sfc_alb_dir_vis`, 0.067, which it used before, left 1003.10 W/m².) At
step 0 Noah-MP has not run yet and its fields hold the undefined placeholder, so the sweep
uses `erf.radiation.surface_albedo_sw`. The numbers
are identical on 1 and 4 ranks. The two-stream coupling to Noah-MP is single-level:
`amr.max_level > 0` stops at start-up.
