# ERF inputs-file generator

`Exec/Tools/erfcreateinput.py` writes a complete, commented ERF inputs file
for a chosen configuration, the dummy auxiliary files the deck names, and a
list of the files that only the user can supply.

## Purpose

A new case needs keys from many parts of ERF, and several defaults stop a run
or leave it without output: `erf.init_type` has no default, `erf.use_gravity`
defaults to false, no plotfile is written without `erf.plot_vars_1`, the log
files stay empty with `erf.v = 0`. The generator writes these explicitly.

* The chosen option of each group is an active line with its sub-keys.
* The other options of the group are commented lines (`#key = value`) beside
  it, each with the keys it needs.
* Values the tool cannot know are marked `USER:` and listed in the summary.
* A combination ERF refuses at start-up is refused with a message naming the
  options; no option is changed silently.

The generator complements the inputs checker, `Exec/Tools/check_erf_inputs.py`. The generator writes a new case;
the checker reads any case, written by hand or by the generator, and reports what would stop the run or is likely
not what was meant, including the contents of the files the deck names (soundings, terrain and roughness tables,
`wrfinput`, `wrfbdy`, `met_em`). Use them in that order: generate, replace the dummy files, then check. The
generator's own tests run every deck it writes through the checker.

The checker reads what the generator writes. On a generated case it reports each dummy file that still holds the
placeholder values (`FILES_TO_REPLACE.txt` carries a sha256 fingerprint of each dummy), each file to supply that is
missing, and the values marked `USER:`. Once a dummy is replaced, its warning goes.

## Files

| File | Role |
|------|------|
| `../erfcreateinput.py` | command line |
| `options.py` | options, choices, the per-setup defaults, presets and case files |
| `presets.py` | the named presets |
| `constraints.py` | option ranges and refused combinations, each with its reason |
| `recipes/grid.py` | problem name, domain, refinement, boundary faces, time step, equation set, advection |
| `recipes/terrain.py` | terrain type and source, vertical stretching, immersed-forcing wall keys, buildings |
| `recipes/physics.py` | turbulence closure, surface layer, forcing, damping and sponge, initial state (idealised and real-data) |
| `recipes/models.py` | moisture, radiation, land surface, forest canopy, wind farm |
| `recipes/output.py` | logs, plotfiles, checkpoints, samplers, post-processing notes |
| `auxfiles.py` | dummy files, each in the layout its ERF reader expects, and the files that cannot be generated |
| `deckwriter.py` | sections and lines of the deck; refuses a key set twice |
| `enums.py` | enum values read from `Source/DataStructs` (`AMREX_ENUM`) |
| `generate.py` | assembly, writing, printed summary |
| `selftest.py` | built-in tests |

## How to run

```bash
python3 Exec/Tools/erfcreateinput.py --list
python3 Exec/Tools/erfcreateinput.py --preset abl_convective_les --outdir my_cbl
python3 Exec/Tools/erfcreateinput.py --preset terrain_site_les --nx 128 --ny 128 --outdir my_site
python3 Exec/Tools/erfcreateinput.py --case-file my_cbl/case.json --heat-flux 0.2 --outdir my_cbl2
python3 Exec/Tools/erfcreateinput.py --setup abl --outdir my_abl
python3 Exec/Tools/erfcreateinput.py --setup abl --turbulence mynn25 --roughness file \
    --surface-heat flux --heat-flux 0.1 --coriolis --forcing geostrophic \
    --output plt,plt2d,chk,profiles,stations --outdir my_pbl
python3 Exec/Tools/erfcreateinput.py --setup abl --terrain file --terrain-method immersed \
    --buildings file --amr 1 --outdir my_site
python3 Exec/Tools/erfcreateinput.py --setup hill --hill schar --dz0 20 --moisture kessler --outdir my_hill
python3 Exec/Tools/erfcreateinput.py --setup wrfinput --moisture morrison --turbulence mrf \
    --radiation rrtmgp --lsm noahmp --outdir my_real
python3 Exec/Tools/check_erf_inputs.py my_abl/inputs
```

A case can be given three ways, which combine: a preset (`--preset NAME`, 22 names listed by `--list`), a case
file (`--case-file FILE`, a JSON object of option names and values), and single options. The case file replaces
the preset's values, and options on the command line replace both. Every case directory gets a `case.json` with
its options, so `--case-file my_case/case.json` writes the same deck again.

| Option | Default | Effect |
|--------|---------|--------|
| `--preset` | unset | a named set of options, for example `abl_neutral_les`, `hill_immersed`, `urban_les`, `real_wrf` |
| `--case-file` | unset | a JSON file of options |
| `--setup` | `abl` | idealised: `abl`, `hill`, `bubble`, `density_current`; real-data: `wrfinput`, `metgrid`. Sets `erf.prob_name`, the domain and the defaults of the options below |
| `--case` | from the setup | `ideal` or `real`; must agree with `--setup` |
| `--init` | from the setup | idealised only: `sounding` (`Input_Sounding`), `uniform`, `isentropic` |
| `--dycore` | `compressible` | `compressible` or `anelastic` |
| `--terrain` | from the setup | `none`, `analytic` (with `--hill woa`, `schar`, `cos4`, `--hill-height`, `--hill-width`), `file`, `file_usgs` |
| `--terrain-method` | `fitted` | `fitted` (`StaticFittedMesh`), `immersed` (`ImmersedForcing`), `eb` (`EB`; ridges with the anelastic set only) |
| `--buildings` | `none` | `file`: buildings by immersed forcing from a table |
| `--dz0` | unset | first-cell height [m] of a stretched vertical grid; positive, below the uniform cell height, and thick enough for ERF's ratio solve to converge. The deck states the ratio ERF will solve |
| `--turbulence` | from the setup | `none`, `dns`, `smagorinsky`, `smagorinsky2d`, `deardorff`, `keqn`, `myj`, `mynn25`, `mynnedmf`, `ysu`, `ysunew`, `mrf` |
| `--surface` | from the setup | `slip`, `noslip`, `most` (`zlo.type = surface_layer`; with immersed terrain, the wall law on the immersed surface) |
| `--roughness` | `const` | `const` (`erf.most.z0`) or `file` (`erf.most.roughness_file_name`) |
| `--surface-heat` | `adiabatic` | `adiabatic`, `flux`, `temp`, or `file`: a sea-surface temperature table over water (`erf.is_land = 0`) |
| `--heat-flux` | 0.1 | kinematic heat flux [K m/s] for `--surface-heat flux` |
| `--moisture` | from the setup | `none`, `kessler`, `kessler_norain`, `satadj`, `moistnocondensation`, `sam`, `sam_noice`, `sam_noprecip_noice`, `morrison`, `morrison_noice`, `wsm6`, `wdm6` |
| `--radiation` | `none` | `twostream`, `rrtmgp` |
| `--lsm` | `none` | `slm`, `noahmp` |
| `--canopy` | `none` | `file`: forest drag from a table of tree patches |
| `--windfarm` | `none` | `fitch`, `ewp`, `simplead` |
| `--amr` | 0 | refined levels, 0 to 2 (x and y refined, ratio 2; y is not refined when `--ny 1`) |
| `--refine` | `box` | with `--amr`: `box`, nested static boxes at the domain centre; or `field`, cells tagged by field tests (erf-model/ERF#5036) |
| `--refine-on` | from the setup | with `--refine field`: one or more conditions `FIELD:TEST:VALUE`. Fields: density, theta, scalar, qv, qc, qi, qr, qs, qg, qt, vorticity, helicity, max_reflectivity, u, v, w, omega, terrain_blanking. Tests: `greater`, `less`, `gradient` (`value_greater`, `value_less`, `adjacent_difference_greater`). Several conditions go under `erf.<indicator>.conditions` |
| `--refine-combine` | `and` | with several conditions: tag where all hold, or `or` where any holds |
| `--refine-box` | unset | with `--refine field`: restrict the tests to `XLO YLO XHI YHI` or `XLO YLO ZLO XHI YHI ZHI` [m] |
| `--refine-z` | `above_ground` (`physical` with immersed or EB terrain) | how the z extent of `--refine-box` is measured (`erf.<indicator>.box_z_coordinate`) |
| `--lateral` | from the setup | `periodic`, `slipwall`, `outflow`, `inflow_outflow`, `open` (the last three act on x; y is periodic); `real` for the real-data setups |
| `--forcing` | `none` | `pressure_gradient`, `geostrophic`, `geo_table` |
| `--coriolis`, `--latitude` | off, 45 | Coriolis force (always on for real-data setups); latitude in [-90, 90] degrees |
| `--damping` | from the setup | `none`, `rayleigh`, `sponge` |
| `--nx --ny --nz` | from the setup | cells; at least 1 in x and y, at least 2 in z. With fewer than 16 cells in y, hills are written as ridges across x |
| `--domain LX LY LZ` | from the setup | extent [m], each positive |
| `--dt` | unset | fixed time step [s], positive; unset gives `erf.cfl` (a fixed step for setups that start at rest) |
| `--stop-time` | from the setup | run length [s], positive |
| `--wind U V` | from the setup | initial and driving wind [m/s] |
| `--theta` | 300 | reference potential temperature [K], positive |
| `--z0` | 0.1 | roughness length [m], positive and below the first cell centre |
| `--output` | `plt,chk,logs` | active outputs among `plt`, `plt2d`, `netcdf`, `chk`, `logs`, `profiles`, `line`, `plane`, `stations`, `subvol`; the others are written as comments |
| `--outdir`, `--name` | `erf_case`, `inputs` | where the case is written; the name must not be one of the generated file names |
| `--force` | off | overwrite an existing case |

An option outside its range, a non-finite number, or an unknown output name stops the script with a message
naming the option. Nothing is written when a case is refused.

Choices the tool makes for the user, each stated in the deck:

* Anelastic decks use the FFT Poisson solve (`erf.use_fft = true`), which needs an FFT build; multigrid is the
  commented alternative.
* The explicit Rayleigh and sponge rates are capped so that rate x time step stays at or below 0.25, and an
  adaptive anelastic deck with damping gets `erf.dt_max`.
* Setups that start at rest get a fixed time step.
* An anelastic bubble perturbs theta, not density (`prob.perturb_rho = false`).
* Field refinement with a turbulence closure sets `erf.implicit_thermal_diffusion` and
  `erf.implicit_momentum_diffusion` to false: tagged boxes can be stacked in z, and ERF stops when the vertical
  solve is implicit on such boxes.

## Dummy files

Each dummy is valid and sized to the grid, so the run starts; its values are
placeholders. `FILES_TO_REPLACE.txt` in the case directory repeats what the
script prints about each one.

| File | Key | Layout | Dummy content |
|------|-----|--------|---------------|
| `input_sounding` | `erf.input_sounding_file` | line 1 `p_surf[hPa] theta[K] qv[g/kg]`; then `z theta qv u v`, z increasing to the domain top | neutral, dry, uniform wind |
| `roughness.txt` | `erf.most.roughness_file_name` | `x y z0` [m], one row per level-0 grid node, x fastest, `(nx+1)*(ny+1)` rows | uniform `--z0` |
| `sfc.txt` | `erf.most.sfc_file` | header line, then `day sst[K]` | constant sea-surface temperature |
| `geo_wind_table.txt` | `erf.abl_geo_wind_table` | `z Ug Vg`, z spanning the domain | uniform wind |
| `terrain.txt` | `erf.terrain_file_name` | one number per line: nx, ny, the x values, the y values, then heights with y fastest | one smooth hill (a ridge with fewer than 16 cells in y), zero at the lateral faces |
| `terrain_usgs.txt` | `erf.terrain_file_name_USGS` | `lon lat`; `nx ny`; then `x y height`, x fastest | the same hill |
| `buildings.txt` | `erf.buildings_file_name` | as `terrain.txt`; heights above the ground | one block at the centre |
| `forest.txt` | `erf.forest_file` | `type xc yc height diameter cd lai laimax` per patch | one patch at the centre |
| `windturbines_loc_x_y.txt` | `erf.windfarm_loc_table` | `x y` [m] per turbine | one turbine at the centre |
| `windturbines_spec.tbl` | `erf.windfarm_spec_table` | row count; `hub_height rotor_diameter standing_Ct nominal_power`; then `wind_speed Ct power` | a 4-row curve, not a real turbine |

## Files the tool cannot generate

These are binary or model data. The deck names them, the summary lists them
under "Files you must supply", and the run stops until they are present.

| Option | Files |
|--------|-------|
| `--setup wrfinput` | `wrfinput_d01`, `wrfbdy_d01` |
| `--setup metgrid` | two or more `met_em` files |
| `--radiation rrtmgp` | the four RRTMGP coefficient files |
| `--lsm noahmp` | `namelist.erf`, `NoahmpTable.TBL` |

For a real-data case `amr.n_cell` and the domain extent must be edited to the
grid of the NetCDF file.

## Refused combinations

Each rule in `constraints.py` carries the source line or the test run behind it.

| Combination | Reason |
|-------------|--------|
| `--surface most` with `--turbulence none` | ERF stops: a surface-layer face needs a diffusive closure |
| a PBL scheme other than `myj` without `--surface most` | ERF stops: the scheme reads the surface-layer fluxes |
| `--turbulence ysu` without a negative `--heat-flux` | ERF stops unless the surface layer is stable |
| `--turbulence keqn` with `--surface slip` | ERF stops: the wall distance needs a solid boundary |
| `--z0` at or above the first cell centre | ERF stops |
| `--dz0` too thin for the cell count | ERF's solve for the stretching ratio does not converge |
| `--forcing geostrophic` or `geo_table` without `--coriolis` | ERF reads the geostrophic wind only with the Coriolis force on |
| `--terrain-method eb` with `deardorff` or `keqn`, `--roughness file`, `--surface-heat file` or `--buildings file` | ERF stops |
| `--canopy file` with `--terrain-method immersed` | ERF stops |
| a `sam` moisture model on terrain-fitted terrain or in a real-data case | ERF stops on a variable-dz mesh |
| `--lsm` without a moisture model | ERF stops |
| `--setup wrfinput` with `--moisture none` | ERF stops |
| `--output profiles` with `--turbulence none` | ERF ended in a segmentation fault at the first profile write |
| `--output plane` with `--amr` | plane sampling over levels needs refinement in z as well |
| `--terrain-method eb` when compressible, over a round hill, or with a terrain table in `--setup abl` | in test runs these stopped at start-up or diverged within 100 steps; an anelastic ridge stayed steady |
| the real-data setups with `--init`, a lateral type other than `real`, terrain options, `--dz0`, `--buildings` or `--forcing` | the NetCDF files set these |
| `--dz0` with `--terrain-method immersed` or `eb`, or with `--buildings file` | those methods keep uniform cells |
| `--lsm` without `--surface most` on the zlo face | ERF prints that the land model will not be applied |
| `--roughness file` or `--surface-heat` other than `adiabatic` without `--surface most`; `--surface-heat file` with `--roughness file` | the keys belong to the surface layer; over water the roughness follows the sea state |
| `--surface-heat` other than `adiabatic`, or `--setup bubble` / `density_current`, with `--init uniform` | a uniform base state runs without gravity, so there is no buoyancy |
| `--lateral slipwall` with a wind; `--lateral inflow_outflow` without a positive x wind; `--setup density_current` with a wind | the wind cannot enter, or the problem does not read it |
| `--hill schar` with fewer than 8 cells per ripple; `--hill cos4` on a ridge grid or wider than the domain | the shape would not be resolved |
| `--terrain-method immersed` with a PBL scheme, `--roughness file` or `--surface-heat file`; `--terrain-method eb` with `--canopy` or `--amr`; `--damping sponge` with a real-data setup | not covered by this tool |
| `--amr` with `--nx` or `--ny` not divisible by 4 (with `--refine box`) | the refined box must lie on coarse cells |
| `--refine field` with a qv, qc, qt, qr, ice or reflectivity field the moisture model does not carry, or terrain_blanking without immersed terrain | ERF stops (`ERF_Tagging.cpp`) |
| `--refine-z above_ground` with a 3-value `--refine-box` over immersed or EB terrain | ERF stops |
| `--refine field` with a PBL scheme | the PBL kernels need boxes over the full column, and tagged boxes need not reach the ground; ERF stopped on an assertion in test runs |
| `--refine field` with `--dycore anelastic` in 3-D | in test runs half of these cases stopped within 20 steps (the pressure solve, an FFT boundary assertion), also after erf-model/ERF#5075; the 2-D density current ran |
| `--refine field` with `--setup bubble` | the Bubble problem needs boxes over the full depth |
| `--refine`, `--refine-on` or `--refine-box` without `--amr`; `--refine-on` or `--refine-box` with `--refine box` | the options belong to field refinement |
| `--output netcdf` without `plt` | `netcdf` is the file type of the plotfile stream |

## Expected results

`python3 Exec/Tools/erfcreateinput.py --self-test` exits 0 with
`erfcreateinput self-test: 157 cases written, 100 refused, 0 run in ERF, 0 failure(s)`;
CTest runs it as `InputsCreator_SelfTest`. The cases include every preset. For each case it checks that

* every key, active or commented, is read under `Source/`, and is not a name the docs list only under another
  prefix;
* the inputs checker reports each dummy file as a placeholder, and no other ERROR or WARNING than the missing
  files the tool cannot generate;
* each dummy file has its reader's layout, with the values of `--theta`, `--wind` and `--z0`;
* the deck parses back to the options, and `case.json` gives the same deck again.

It also checks that each refused combination is refused with a message naming the option, and that a
command-line option replaces a preset's or case file's value.

Optional, slower checks:

* `--self-test --schema schema.json` also requires every key to be in a schema written by
  `build_erf_schema.py`. Known exceptions, for the reason given in erf-model/ERF#5067 (keys read through a name
  held in a variable): four SLM soil keys, `value_less` and `adjacent_difference_greater` of a refinement
  indicator, and the keys of a refinement condition (`erf.<indicator>.<condition>.*`).
* With `ERF_CREATEINPUT_EXE=/path/to/erf_exec` each case whose files are all generated is run for three steps
  with every active output written each step; the run must finish and AMReX must list no unused input.
  `ERF_CREATEINPUT_FILTER=text` runs only the cases whose options hold that text, and
  `ERF_CREATEINPUT_NETCDF=1` includes the NetCDF plotfile case.

Runs made on 2026-11-12 with ERF at development 80427a4c3 (Release, MPI, FFT, NetCDF), which includes the field
tagging of erf-model/ERF#5036 and the fixes #5070, #5072 and #5075:

| Build | Cases | Result |
|-------|-------|--------|
| no wind farm | 145 (all but the NetCDF plotfile case and the cases waiting for supplied files) | all ran |
| wind farm on | 4 wind-farm cases, and a case without a wind farm | all ran |
| NetCDF and wind farm on, at abce396eb | a NetCDF plotfile case; a Metgrid case on `Exec/RegTests/MetGrid/met_em.d01.*`; the WRF-input deck on `Exec/CanonicalTests/EkmanSpiral/wrfinput_ekman_d01`, with periodic sides in place of the boundary file | all ran |

Longer runs of single cases (30 to 100 steps) were made for the anelastic bubble, the sponge and Rayleigh
rates, the isentropic boundary layer and the EB ridge; the field-refinement cases ran 12 to 20 steps with a regrid
every 2 steps. The three review passes ran about 430 further cases.

Not run: a WRF-input case with a wrfbdy file (none is in the repository), RRTMGP and Noah-MP (their builds and
data were not available). Those decks are checked against the schema and the inputs checker only.

ERF behaviours found with the generated cases, and where they stand:

| Behaviour | Issue | Now |
|---|---|---|
| anelastic multigrid did not converge on cells with aspect ratios of 4:1 and more | erf-model/ERF#5064 | fixed by erf-model/ERF#5075 (coarsening per direction): the three failing cases ran with multigrid at 80427a4c3; the decks still use the FFT solve |
| a build with `ERF_ENABLE_WINDFARM` stopped on every case without a wind farm | erf-model/ERF#5066 | fixed by erf-model/ERF#5070 |
| EB with `zlo.type = surface_layer` ended in a segmentation fault unless `erf.eb_boundary_type = SurfaceLayer` | erf-model/ERF#5068 | ERF stops at start-up with a message since erf-model/ERF#5072; the decks set the type |
| five SLM keys and the keys read through a name in a variable are missing from the inputs schema | erf-model/ERF#5067 | open; listed as exceptions in the schema check |
| the profile logs (`erf.profile_int`) ended in a segmentation fault at the first write after step 0 without a turbulence closure | not filed | refused by the tool |
| EB terrain stopped at start-up when compressible on the default hill grid, and diverged within 100 steps over a round hill and for a table ridge in the boundary-layer setup | not filed | refused by the tool except for an anelastic ridge |
| field tagging: the PBL kernels assert that every box spans the column; half of the 3-D anelastic cases stopped within 20 steps | not filed | refused by the tool |

## Adding an option

1. Add the choice to `options.py` (and a default to `SETUP_DEFAULTS` when it
   depends on the setup).
2. Write its keys in the matching recipe with `Section.set`, `Section.user`,
   `Section.alt` or `Section.choice`. Keep the other choices as alternatives.
3. Add a rule to `constraints.py` for every combination ERF refuses, citing
   the source line.
4. If it names a file, add a writer to `auxfiles.py` and a layout check to
   `selftest.py`; a file that cannot be generated goes in `required_files`.
5. Add the case to `CASES`, and each refused combination to `REFUSED`, in
   `selftest.py`, and run the self-test with `ERF_CREATEINPUT_EXE` set.

## Not included

SuperDroplets and spectral-bin microphysics, SHOC, GeneralAD wind farms, STL
buildings, NetCDF terrain and gridded forest files, moving terrain, boundary
plane input and output, nudging and large-scale forcing files, hindcast
initialisation, particles, the fire and moving-body couplings.
