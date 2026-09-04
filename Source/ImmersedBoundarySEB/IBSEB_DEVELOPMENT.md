# Immersed-boundary surface energy balance (IBSEB): development plan

Branch `claude-seb`, based on `development`. Nothing is merged in from the
other branches of this fork; the pieces that are reused (`ERF-SLUCM`
kernels, `ERF-Radiation` solar geometry) are copied as single files, so the
branch carries only its own commits and the PR to `erf-model/ERF`
`development` is reviewable on its own whatever the state of those
branches. Radiation reaches the balance through a small provider interface
(below), so the branch builds and tests without the two-stream model and
gains it when that model is available.

The balance runs on the faces of resolved buildings (immersed forcing,
`erf.buildings_type = ImmersedForcing`), one skin temperature and one slab
per face:

    C dT/dt  =  SW_abs + LW_net - H - LE - G          (+ Q_ext, an external
                                                        incident flux hook)

Every term is either reused from code that exists on a branch of this fork
or replaced by resolved geometry:

| term | source |
| --- | --- |
| skin-temperature Newton solve, slab conduction, material library | `ERF-SLUCM` kernels (`ERF_UCMSEBSolver.H`, `ERF_UCMSlabConduction.H`, `ERF_UCMMaterialRegistry`), lifted verbatim with the canyon arguments renamed |
| downwelling shortwave and longwave at the face height | a radiative-forcing provider: prescribed (this branch) or the `ERF-Radiation` two-stream column (when that branch is merged) |
| sun position | `ERF-Radiation` `ERF_SolarGeometry.H`, copied as one file |
| sensible heat through the wall | the immersed-forcing MOST wall model already in `ERF_ImmersedForcing.cpp`, generalised from the top face to all three face directions and driven by the face skin temperature |
| shadowing, sky / ground / building view fractions | new: ray casting against the building height map |
| face storage | new: compact per-rank face list, owned by the adjacent fluid cell |

What is *not* built: canyon morphology, Facet3D injection, BEP drag, the
Kusaka canyon-air budget, canyon view factors, radiosity. The immersed
boundary resolves all of that.

## Radiative forcing provider

The balance asks one interface for what it needs at a face:
`SW_direct_down(z)`, `SW_diffuse_down(z)`, `LW_down(z)` [W/m2] and the unit
sun vector. Two providers:

- **prescribed** (`erf.ibseb.radiation = prescribed`): sun vector from the
  copied solar geometry at the site and time, or a fixed zenith and azimuth
  for tests; clear-sky direct and diffuse from the Bird transmission form
  and a gray-sky longwave from the air temperature, both height-independent.
  This is what the regtests of phases 2, 3 and 6 use, and what a run without
  a radiation model gets.
- **two-stream** (`erf.ibseb.radiation = two_stream`): the column fluxes of
  `ERF-Radiation` at the face's own level, with the direct beam stored as
  a fifth component of `rad_fluxes` (the column keeps it internally but
  stores only the sum today). Built on a local merge of `ERF-Radiation`
  into a working copy when Phase 6 wants a real diurnal cycle, kept behind
  the interface, and included in the PR only if that branch has merged
  upstream by then; otherwise it follows as a small second PR.

## Conventions

- Opt-in: `erf.ibseb.enable = false` by default; nothing allocated and no
  result changes without it.
- Every phase ends with a regtest under `Exec/RegTests/IBSEB_<phase>` that
  runs on one and four ranks and gives identical numbers, plus a README with
  the reference values, and a section in this file.
- Signs follow the SLUCM / WRF convention: positive H, LE, G leave the
  surface; positive SW_abs and LW_net enter it.
- One update per atmospheric step, before the MRI integrator, like the
  land-surface models; the face flux is added (`+=`) to the cell source
  after `make_sources`, never overwritten.
- GPU: all per-face work is `ParallelFor` over the face list; host loops
  only at initialisation and for output.

## Storage (Phase 1)

`IBFaceSet` per level, built once after `init_immersed_forcing`:

- A face is a wall face where the face blanking is 1 on the solid side and
  0 on the fluid side (x, y or z faces; z faces with fluid above are roofs,
  with fluid below are not walls). Faces whose fluid cell lies in this
  rank's valid region belong to this rank, so every input to the balance is
  local and the per-step update needs no communication.
- Struct of arrays in `amrex::Gpu::DeviceVector`: `i, j, k` of the fluid
  cell; `dir` (0/1/2) and `side` (which side the solid is on); `bid`
  (building id, from the height map's connected footprints, the labelling
  the fire module already has); `mat` (material id); `area`; `T_skin`;
  `T_slab[n_layers]`; `f_sky, f_ground, f_bldg` (Phase 3); `SW_abs, LW_net,
  H, LE, G, Q_ext` (current fluxes); `n_hits` scratch.
- Replicated on every rank as device vectors: the building height map (top
  height per column) and terrain height, for ray casting. A 1000 x 1000
  column domain is 8 MB per rank.
- Output: scatter into face-centred `MultiFab`s on demand (`ibseb_tskin_x`,
  `_y`, `_z`, and the fluxes) for the plotfile; checkpoint the same way and
  refill the list from the fields on restart, so restarts are independent
  of the rank count. Per-building CSV (mean and max skin temperature by
  face orientation, fluxes) every `erf.ibseb.csv_int` steps.
- Inputs: `erf.ibseb.enable`, `erf.ibseb.n_slab_layers` (4), `erf.ibseb.
  material_file` (CSV in the SLUCM schema), `erf.ibseb.material_default`,
  `erf.ibseb.T_interior` (fixed interior temperature), `erf.ibseb.T_skin_init`,
  `erf.ibseb.csv_file`, `erf.ibseb.csv_int`, `erf.ibseb.plot_faces`.
- Regtest `IBSEB_Storage`: the ImmersedForcingTest skyscraper on 1 and 4
  ranks; face counts per direction and per building, sum of face areas
  against the analytic box area, a plotfile with the scattered skin
  temperature, and a checkpoint / restart round trip that reproduces the
  list bit for bit.

## Shortwave (Phase 2)

- Sun vector and the direct and diffuse beams from the provider (Phase 2
  builds the prescribed one; the two-stream one is Phase 6's option).
- Direct on a face = beam x max(0, n . s), zero when the ray from the face
  centre toward the sun hits a building or the terrain. Ray cast: a 2D DDA
  over columns of the height map with the ray height rising along it.
- Diffuse on a face = f_sky x diffuse down + f_ground x ground albedo x
  (direct + diffuse) at the ground. f_sky and f_ground come from Phase 3's
  hemisphere sampling; in Phase 2 they are placeholders (1 for roofs, 0.5
  for walls) and the regtest is on the direct beam.
- Absorbed = (1 - albedo_mat) x incident.
- Regtest `IBSEB_Shortwave`: one 20 m cube at a prescribed sun; sunlit and
  shaded faces against the analytic incidence angles, the shadow length on
  the ground behind the cube against H / tan(elevation), and a closure
  check that the absorbed sum never exceeds the incident sum.

## Longwave (Phase 3)

- Hemisphere sampling at initialisation: N_az x N_el rays per face (32 by
  default) through the same ray caster; each ray ends in sky, ground or a
  building, giving `f_sky + f_ground + f_bldg = 1`.
- LW_in = f_sky x LW_down at the face height (provider) + f_ground x
  eps_g sigma T_ground^4 + f_bldg x eps sigma T_skin^4, the last being the
  isothermal-surroundings approximation (a face and the walls it sees
  exchange no net longwave). T_ground is the surface layer's surface
  temperature at the column the ray hits, or its mean under the building.
- LW_net = eps x (LW_in - sigma T_skin^4).
- No face-to-face view factors and no radiosity in this plan; the
  fractions are stored so a radiosity pass can be added later without
  touching the balance.
- Regtest `IBSEB_Longwave`: closure of the three fractions to one on every
  face; a night case with fixed skin temperatures where the face-mean
  longwave loss matches the analytic value for an isolated cube.

## Sensible and latent heat (Phase 4)

- Sensible: the immersed-forcing MOST wall model already computes a
  friction velocity and a target temperature from a surface temperature
  and the tangential wind at the first fluid cell. Generalise
  `ImmersedForcingBuildings_Scalar` from "solid below, fluid above" to all
  three face directions with `T_skin` of the face as the surface
  temperature; neutral similarity on vertical faces, the existing stability
  correction on roofs. H per face is diagnosed from the same u*, theta*
  pair so the balance and the atmosphere see one flux.
- The forcing keeps its relaxation form (it is what the anelastic and
  compressible paths both support); the equivalent face flux is stored on
  the face for the balance and for output.
- Latent: `LE = 0` with the argument in place; a wet-surface option
  (`erf.ibseb.wet_fraction` per material, bulk evaporation) is a later
  addition and needs the moisture source slot.
- Regtest `IBSEB_Sensible`: a cube with a prescribed hot skin in a neutral
  channel; face fluxes against the MOST formula evaluated offline from the
  plotfile wind, and the domain heat budget (sum of H x area over faces
  equals the rate of change of integrated rho cp theta) to a few percent.

## Ground heat flux (Phase 5)

- Slab conduction per face with the SLUCM implicit tridiagonal solver:
  `n_slab_layers` layers of material thickness / n each, flux boundary at
  the skin, fixed `T_interior` at the inside. G is the flux into the slab.
- Materials by building id from the material CSV (albedo, emissivity,
  conductivity, heat capacity, thickness); one default material otherwise.
- Interior temperature fixed per run in this plan; a lumped interior
  budget per building (sum of inward conduction over its faces) is the
  follow-up that the building ids make possible.
- Regtest `IBSEB_Ground`: a step change in skin temperature on a thick slab
  against the semi-infinite erfc solution; a thin slab against the lumped
  exponential.

## Prognostic balance (Phase 6)

- Per face and per step: gather SW_abs, LW_in, wind and air temperature at
  the fluid cell, T1 of the slab; Newton on T_skin with the SLUCM solver
  extended by `Q_ext` (an external incident flux, absorbed with the face's
  absorptivity, zero here; the fire's radiant flux later enters through
  it); advance the slab with the resulting G; store the fluxes.
- Longwave emission linearised inside Newton as it already is; the
  20 K step clamp and the 260-380 K bounds stay, with the bounds as
  inputs since fire exposure will exceed 380 K.
- Diagnostics: per-building CSV, `[IBSEB]` summary line (faces, min / max
  skin temperature, largest residual), plotfile faces, checkpoint of
  `T_skin` and the slab.
- Optionally the two-stream provider, on a local merge of `ERF-Radiation`,
  behind the interface.
- Regtest `IBSEB_Prognostic`: 24 h on one cube with the prescribed diurnal
  cycle; the balance residual on every face below a tolerance
  every step, the energy budget closed (radiation in = stored + convected +
  conducted), restart reproduces the straight run, and a comparison of the
  roof against the ground surface layer's own temperature, which sees the
  same sky.

## Canonical: isolated building (Phase 7)

`Exec/CanonicalTests/IBSEB/IsolatedBuilding`: a 20 m cube on a 5 m grid,
neutral sounding, MOST ground, prescribed clear-sky radiation (two-stream
where available), 24 h from sunrise.
Documented: the sunlit / shaded wall temperature contrast through the day
(order 10-20 K for masonry at midday), the roof against the ground, the
phase lag of the wall temperature behind the sun, and the wake-side wall
cooling. `check_isolated_building.py` reads the CSV and asserts the
qualitative sequence and the closure.

## Canonical: building set (Phase 8)

`Exec/CanonicalTests/IBSEB/BuildingSet`: the three-box street from the
WUI obstacle decks and the skyscraper block from the immersed-forcing
tests. Exercises mutual shadowing, the building fraction of the view
fractions, several materials, and the per-building CSV. Documented with
timing (faces per rank, ray-cast cost at initialisation and per step) so
the cost of the face list on a city-scale case can be estimated.

## After the eight phases

- PR to `erf-model/ERF` `development`, self-contained; the two-stream
  provider goes with it or follows it depending on the radiation branch.
- ERF-Hazard: merge development, then the fire's radiant flux into
  `Q_ext`, the exposure metrics from the skin temperature, ignition, and
  the drop of the atmospheric injection to the convective fraction.
- Later options, none needed for the PR: face-to-face radiosity, lumped
  interior budget, wet surfaces, terrain faces under the same balance.

## Phase log

| phase | status | regtest | notes |
| --- | --- | --- | --- |
| 1 storage | planned | `IBSEB_Storage` | |
| 2 shortwave | planned | `IBSEB_Shortwave` | prescribed provider, ray-cast shadow |
| 3 longwave | planned | `IBSEB_Longwave` | |
| 4 sensible, latent | planned | `IBSEB_Sensible` | generalises the immersed-forcing wall model |
| 5 ground | planned | `IBSEB_Ground` | lifts the SLUCM slab solver and materials |
| 6 prognostic | planned | `IBSEB_Prognostic` | lifts the SLUCM Newton solver, adds `Q_ext`; two-stream provider optional |
| 7 isolated building | planned | canonical | |
| 8 building set | planned | canonical | |
