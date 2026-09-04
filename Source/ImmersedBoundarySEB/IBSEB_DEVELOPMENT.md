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
- Every phase ends with a regtest under `Exec/CanonicalTests/SEB/<phase>` that
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
- Regtest `SEB/Phase1_Storage`: the ImmersedForcingTest skyscraper on 1 and 4
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
- Regtest `SEB/Phase2_Shortwave`: one 20 m cube at a prescribed sun; sunlit and
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
- Regtest `SEB/Phase3_Longwave`: closure of the three fractions to one on every
  face; a night case with fixed skin temperatures where the face-mean
  longwave loss matches the analytic value for an isolated cube.

## Sensible and latent heat (Phase 4)

- Sensible: a wall function on every face with the same form as the
  immersed forcing's momentum wall model: u* from the tangential wind of the
  fluid cell at half a cell from the wall with z0_wall, theta* from the
  skin-to-air potential-temperature difference with z0h_wall, H = rho c_p u*
  theta*; neutral on walls, the surface layer's similarity functions on
  roofs when asked.
- The flux enters as an explicit source (H A / (c_p V Pi) into rho-theta of
  the fluid cell, added after make_sources at every slow stage), chosen over
  the relaxation form because it is exact and local; the immersed forcing's
  own surface-temperature inputs are refused when the balance is on.
- Latent: `LE = 0` with the argument in place; a wet-surface option
  (`erf.ibseb.wet_fraction` per material, bulk evaporation) is a later
  addition and needs the moisture source slot.
- Regtest `SEB/Phase4_Sensible`: a cube with a prescribed hot skin in a neutral
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
- Regtest `SEB/Phase5_Ground`: a step change in skin temperature on a thick slab
  against the semi-infinite erfc solution; a thin slab against the lumped
  exponential.

## Prognostic balance (Phase 6)

- Per face and per step, after the three flux routines: the implicit slab's
  linear response `G = a T_s - b` to the skin temperature (two trial slab
  steps, ERF_IBSEBSlab.H `slab_skin_response`), then Newton on
  `SW_abs + eps Q_ext + eps [LW_ext - (1 - f_b) sigma T_s^4] - C_H (T_s/Pi -
  theta_a) - LE - (a T_s - b) = 0` (ERF_IBSEBBalance.H, lifted from the
  SLUCM `solve_facet_seb`), then the slab advanced with the solution and
  LW_in, LW_net, H, G rewritten at it. The skin is massless; the slab's
  first layer carries the storage, as in the SLUCM.
- Why the response and not the lagged `2k (T_s - T_0^n)/dz`: with the lag
  the flux the balance closes with and the flux the slab absorbs differ by
  one step, and the closure check would show it. With the response they
  agree to rounding and the slab energy changes by exactly
  `dt (G - G_bottom)`.
- `Q_ext` is an *incident* flux absorbed with the longwave emissivity
  (Kirchhoff; fire radiation is thermal). `erf.ibseb.Q_ext_uniform` sets it
  on every face for tests; the fire coupling will write `d_Q_ext` directly
  before `solve_balance`.
- The wall function's coefficient `C_H` is frozen at the wind of the step
  (with the stability functions, at the previous skin temperature). The
  wall term of the incoming longwave folds into the emission as `(1 - f_b)`.
- Clamps: `T_skin_min` / `T_skin_max` (260 / 380 K) and `newton_max_step_K`
  (20 K) as inputs; a face at a bound keeps a non-zero residual, stored per
  face (`d_resid`, `d_niter`), in the summary (`resid_max`), the CSV and
  the dump.
- `erf.ibseb.prognostic = false` keeps the fixed skin of phases 2-5; their
  decks set it.
- Diagnostics: `[IBSEB]` line gains `Q_ext_mean` and `resid_max`; the CSV
  gains `G_mean_Wm2, Q_ext_mean_Wm2, T_skin_min_K, T_skin_max_K,
  resid_max_Wm2`; the dump gains `H_coeff, Q_ext, LE, LW_ext, resid, n_iter`
  and every slab layer; `erf.ibseb.dump_faces_tag_step` keeps one dump per
  step. The checkpoint is unchanged (skin and slab).
- Regtest `SEB/Phase6_Prognostic`: a 40 m cube on a 320 m, 10 m grid; the
  closure deck (fixed sun, 200 steps, every step dumped): residual < 1e-3
  W/m2 on every face, stored H / LW_net / G consistent with T_skin to 1e-6,
  slab energy per step exact to the dump precision, closure over the run
  within the summed residual, and an independent Python model (own Newton,
  dense implicit slab) driven by the dumped forcing reproducing T_skin to
  1e-9 K; the external-flux deck (3000 W/m2, bound raised to 700 K, the
  faces past 380 K, same closure checks); checkpoint restart; and the
  sunrise deck (solar mode, Boulder, solstice, 90 min at 1 s) where the
  east wall warms before the roof and the west wall.
- Not done here: the two-stream provider (waits on ERF-Radiation merging;
  the interface is in place), the roof-versus-ground comparison (the
  development branch's land-surface models have no radiation to see the
  same sky with), and a full diurnal cycle, which is the Phase 7 canonical.

## Canonical: isolated building (Phase 7)

`Exec/CanonicalTests/SEB/Phase7_IsolatedBuilding`: a 40 m cube as an exact
box (`eb2.geometry = box`, faces on cell faces, so 80 plane faces and no
partial cells) on the phase 6 grid (10 m, 320 m periodic domain so a day
runs in about 100 minutes on four ranks), 30 cm of concrete in ten layers, a MOST ground, a 3 m/s
westerly, neutral at 300 K; Boulder on the June solstice from midnight
solar time for 24 h with the prescribed clear-sky provider and a gray sky.
The per-building CSV gains the sun (zenith, azimuth, DNI, diffuse) so the
day can be plotted without debug output; the faces are dumped every
minute into `faces/`.

`check_isolated.py` asserts: the balance residual over the day; night
cooling of every face below the air with a negative net longwave; the
east wall warming first after sunrise; the peaks in the order east, roof,
west with the roof in the early afternoon; a south-north contrast at
midday; the daily absorbed shortwave on the core roof against the
clear-sky formulas integrated independently in Python; and the slab
energy over the day against the integrated conduction from the dumps.
`plot_isolated.py` draws the skin temperature by orientation with the air,
the roof budget, the sun path, a slab Hovmoller of the roof and yt slices.

Not in this case: the ground's own balance (the development branch's
land-surface models have no radiation, so the ground stays at 300 K for
the walls' longwave and MOST drag only), and the two-stream provider.
A 1 s step traps the dycore's floating-point check on this grid with the
balance off as well (noted under findings), so the day runs at 0.5 s.

## Canonical: building set (Phase 8)

`Exec/CanonicalTests/SEB/Phase8_BuildingSet`: the three-box street from the
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

## Findings outside the balance

- The immersed-forcing atmosphere of `development` does not restart bit-for-bit
  (Phase5_Ground, thin deck: the wind at the faces differs by about 1e-4 after
  a checkpoint restart; the slab itself restarts exactly). To be reported
  separately once isolated.
- Phase 6 caught a sign error of phase 2: `solar_azimuth` had the sine
  term with the wrong sign, mirroring the sun east-west (morning sun in the
  north-west). The noon check of phase 2 could not see it (azimuth 181
  instead of 179); the sunrise deck's "east wall warms first" did. Fixed in
  phase 6; the phase 2 check now expects the azimuth just under 180 at
  12:00 MST.
- `erf.fixed_dt = 1.0` on the 10 m, 3 m/s sunrise deck trips the invalid
  floating-point trap in the dycore at step 2 (with the balance converged
  and the skin at 299.5-299.8 K); 0.5 s runs fine. Not investigated; the
  deck uses 0.5 s.
- The immersed forcing goes unstable on the partial cells of a nodal
  height-map building: the phase 6 cube (40 m on a 10 m grid) has corner
  cells with 1-20 % solid from the reader's one-cell ramp, and in the 24 h
  deck a vertical 2-cell checkerboard in theta (+-10 K, w of 5 m/s) grows
  in one such column and trips the floating-point trap after 1.8-2.5 h,
  with the balance off, with a no-slip ground, with the forcing outside
  the substeps (negative density at 3 min) and with `eb2.small_volfrac =
  0.25` (trap at step 14); halving the step to 0.25 s survived the same
  2.6 h window, so the growth depends on the step and is not a plain CFL
  violation. Reading the buildings forcing: the drag branch scales with
  the fraction and vanishes as it goes to zero, but the MOST branch (a
  partial cell whose normal neighbour is exactly fluid) relaxes at a rate
  independent of the fraction toward `(1 - f)` times the log-law velocity,
  and the temperature branch relaxes toward the "inside" neighbour at a
  fixed rate; so the top cell of a 1 % stub is forced like a full wall on
  top of an almost free cell, a one-cell shear layer and a 2-cell mode.
  An exact box (`eb2.geometry = box` on cell faces) has no partial cells;
  the canonical uses it. Fixed on this branch behind
  `erf.if_snap_partial_cells` (default false, the original selection): the
  wall law and the surface conditions sit on cells at least half solid
  whose normal neighbour is less than half solid, the same rule the balance
  uses for its faces; regtest `Exec/RegTests/ImmersedForcingTest/PartialCells`.

## Phase log

| phase | status | regtest | notes |
| --- | --- | --- | --- |
| 1 storage | done 2026-09-03 | `SEB/Phase1_Storage` | 2056 faces on the skyscraper, rank-independent, restart round trip |
| 2 shortwave | done 2026-09-03 | `SEB/Phase2_Shortwave` | prescribed provider, ray-cast shadow; flag matches an independent cast on all 2616 faces |
| 3 longwave | done 2026-09-03 | `SEB/Phase3_Longwave` | hemisphere sampling, sky/ground/building terms; fractions match an independent sampling on every face |
| 4 sensible, latent | done 2026-09-03 | `SEB/Phase4_Sensible` | wall function per face, explicit face flux into the rho-theta source; internal-energy budget closes against the diagnostic run |
| 5 ground | done 2026-09-03 | `SEB/Phase5_Ground` | slab per face (SLUCM solver with a skin-temperature top), materials by building; erfc and steady checks exact |
| 6 prognostic | done 2026-09-04 | `SEB/Phase6_Prognostic` | Newton balance consistent with the implicit slab, `Q_ext` hook, clamps as inputs; closure and an independent re-integration to 1e-9 K |
| 7 isolated building | in progress 2026-09-04 | `SEB/Phase7_IsolatedBuilding` | 24 h at Boulder on the solstice; sequence, residual, roof shortwave and slab energy checks |
| 8 building set | planned | canonical | |
