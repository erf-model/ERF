# Topographic shading with RRTMGP and Noah-MP

This variant of `NoahMP_Ideal` places a radial Witch of Agnesi hill (250 m high,
250 m half-width) in a 6.4 × 3.2 km domain with 100 m cells. A fixed sun 10° above
geographic east casts a westward shadow. The atmosphere extends to 16 km on a
static terrain-following mesh. Both horizontal directions are periodic.

Build ERF with `ERF_ENABLE_MPI=ON`, `ERF_ENABLE_NETCDF=ON`,
`ERF_ENABLE_NOAHMP=ON`, and `ERF_ENABLE_RRTMGP=ON`. The NetCDF installation must
include parallel I/O and Fortran support. Preparation/checking requires Python,
NumPy and netCDF4; plotting also requires Matplotlib.

From the repository root:

```sh
python3 Exec/RegTests/NoahMP_Ideal/TopographicShading/prepare.py /tmp/shading-hill
cd /tmp/shading-hill
mpiexec -n 2 /path/to/erf_exec inputs > simulation.log 2>&1
python3 /path/to/ERF/Exec/RegTests/NoahMP_Ideal/TopographicShading/check.py \
    surface00001_d01.nc --held surface00002_d01.nc --plot shadow.png
```

`prepare.py` requires an empty run directory. It enlarges the existing ideal
land fixture, assigns the hill's surface height to HGT, sets the Noah-MP timestep
to 1 s, and links the four radiation tables from the RRTMGP submodule. It does
not alter the original `NoahMP_Ideal` case. A successful simulation produces
`surface00003_d01.nc` and `plt2d00003`, not just an exit status: Noah-MP's Fortran
checks can terminate a failing run before those outputs exist.

The native `plt2d*` files can be opened in ParaView or VisIt. Select
`topographic_shadow` and overlay `z_surf` contours. The checker uses the terrain
in the NetCDF output to calculate an independent eastward horizon and compares
every mask cell. It also checks that Noah-MP received `terrain_sw_flux_dn`, that
diffuse sunlight remains in shadow, that sun-facing slopes receive enhanced
shortwave, and (with `--held`) that a skipped radiation update does not apply the
correction again. The final output can be checked with the same command, without
`--held`.

With `ERF_ENABLE_TESTS=ON` the case is also registered as
`ctest -R '^NoahMP_TopographicShading$' --output-on-failure`, when Python
NumPy/netCDF4 are available. It needs no stored gold plotfile.

## Controls and outputs

Both switches default to false:

- `erf.topo_shading = true`: remove direct sunlight obstructed by resolved terrain.
- `erf.slope_rad = true`: project direct sunlight onto the terrain normal.
- `erf.shadlen = 25000`: maximum horizontal search distance in metres.
- `erf.fixed_solar_azimuth_angle`: degrees clockwise from geographic north,
  required when topographic radiation uses a positive
  `erf.fixed_solar_zenith_angle` (the existing parameter is a **cosine**, not degrees).
  Omit both fixed-sun parameters for the astronomical sun.

Available 2D outputs:

| Field | Meaning |
| --- | --- |
| `topographic_shadow` | 1 when other terrain obstructs direct sunlight; 0 otherwise, including night |
| `terrain_direct_factor` | Multiplier applied to horizontal direct-beam flux |
| `terrain_incidence` | Nonnegative cosine between the sun and the terrain normal |
| `unshaded_sw_flux_dn` | Original downwelling surface shortwave from RRTMGP |
| `terrain_sw_flux_dn` | Terrain-adjusted downwelling shortwave supplied to land models |

Before the first radiation update these diagnostics are `-999`. A back-facing
slope can have zero direct flux without a cast shadow; inspect
`terrain_direct_factor` as well as the mask. To visualize the mask by itself,
use `erf.slope_rad=false` and keep shading enabled.

## Scope

The correction currently supports RRTMGP on `StaticFittedMesh`, including grids
initialized through WRFInput and Metgrid. It takes terrain from the finalized
bottom nodal coordinates and uses the grid's map factors and latitude/longitude
rotation. No additional WPS terrain fields are required. It modifies surface
shortwave, preserving diffuse and longwave radiation and the original solar
zenith cosine. RRTMGP's atmospheric column fluxes and heating remain unchanged,
as in WRF's surface terrain correction. Spectral direct/diffuse fields are
handled individually for SLM; Noah-MP receives their broadband sum. TwoStream
is not yet connected and enabling these options with it is rejected.

Terrain patches include off-rank data, periodic images, and ancestor terrain
outside fine AMR grids. Raw radiation is interpolated before applying each
level's terrain correction. Terrain caches are rebuilt after regridding; the
first radiation call after restart or grid creation refreshes missing raw
forcing even between scheduled updates. Horizontal boxes must be disjoint
(`amr.no_box_split_dir=2`); moving meshes are not supported. Fine Noah-MP levels
retain their existing rule: a level without its own land setup file inherits
level-0 land state rather than running its own land model.

Rays use a local tangent plane with the target's map factors over the search
radius; outside a nonperiodic domain there is no terrain information and the
ray stops. Coarse terrain outside nested grids is prolonged piecewise constantly.
The resolved terrain is interpolated between mass points along the ray. For a
dynamic sun, the shadow and slope use the midpoint of RRTMGP's radiation
interval; fields are held until the next radiation update. Use short intervals
when examining moving shadows. Near the horizon (cosine ≤ 1e-4) corrections are
bypassed to avoid division by a tiny number. Diffuse sky-view reduction and
terrain-reflected radiation are not included.

## Validation of this implementation

The hill case passed with RRTMGP + Noah-MP on one and two CPU MPI ranks. All
44 shadowed cells match the independent horizon calculation; the mask, radiation
and Noah-MP outputs agree between decompositions at all three steps. A restart
at step 1 also passed, including a forced radiation refresh on the normally
skipped step. The 10 targeted geometry/forcing tests passed on one and two ranks,
and the RRTMGP/Noah-MP build passed all 1,226 tests carrying the `unit` label.
GPU execution and complete WRFInput/Metgrid real-case runs have not been tested
for this change.
