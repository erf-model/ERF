# ERF inputs checker

`Exec/Tools/check_erf_inputs.py` reads an ERF inputs file the way ParmParse
does (comments, `\` continuations, `FILE =` includes resolved from the run
directory, `#if AMREX_SPACEDIM` / `#ifdef AMREX_USE_GPU` branches, Fortran
namelists skipped, several `key = value` per line, quoted strings, last
definition wins) and reports what would stop the run, what is very likely not
what the author meant, and what could be chosen better. Nothing is executed.

## Purpose

ParmParse never complains about a key that nobody reads, and several of the
ERF start-up checks only fire minutes into a real-data run. The checker finds
those problems from the deck alone:

| Level   | Meaning                                                              |
|---------|----------------------------------------------------------------------|
| ERROR   | the run stops at start-up (an AMReX or ERF abort), or a key cannot be read |
| WARNING | the run starts but a setting is almost certainly not what was meant, or is known to go unstable |
| ADVICE  | a working setting that could be chosen better                        |
| INFO    | derived numbers: cell sizes, estimated time step and substeps, box counts |

## Files

| File | Role |
|------|------|
| `../check_erf_inputs.py` | command line: globs, options, exit code, `--self-test` |
| `deck.py` | ParmParse-compatible parser; `Deck` with typed getters, duplicates, file resolution |
| `context.py` | facts shared by the checkers: cell sizes, levels, anelastic or compressible, idealised or real-data, periodicity |
| `findings.py` | `Finding`, `Severity`, `Report` (text and JSON) |
| `checker.py` | `Checker` base class, registry, `run_all_checks` |
| `checks/keys.py` | duplicate keys, retired keys and values, misspelt or wrongly cased keys (ParmParse keys are case-sensitive; compared with the Sphinx docs and the strings read in `Source/`), `FILE =` found only beside the deck |
| `checks/domain.py` | required keys, cell counts, cell sizes, aspect ratio |
| `checks/grid.py` | blocking factor and `max_grid_size` rules AMReX enforces, odd and uneven boxes, vertical decomposition, boxes per rank |
| `checks/refinement.py` | indicators, boxes, ratios, coupling |
| `checks/timestep.py` | `TimeStepChecker` (shared), `CompressibleTimeStep` (acoustic substeps), `AnelasticTimeStep` (advective limit) |
| `checks/boundary.py` | face types, inflow data, periodic conflicts, surface layer |
| `checks/init.py` | `IdealInitChecker` (soundings, gravity), `RealInitChecker` (WRFInput / Metgrid / NCFile) |
| `checks/physics.py` | turbulence closures, advection schemes, forcing and damping, moisture, terrain and stretching; enum names read from `Source/DataStructs` (`AMREX_ENUM`), exact case where ERF compares strings |
| `checks/output.py` | plotfiles, checkpoints, restarts |
| `checks/files.py` | every file the deck names: exists, and when present its contents against the deck (wrfinput, wrfbdy, wrflow, met_em, terrain, soundings, tables) |
| `selftest.py` | 109 cases, each a deck that must trip a check and a sibling that must not |

## How to run

```bash
python3 Exec/Tools/check_erf_inputs.py Exec/ABL/inputs_most
python3 Exec/Tools/check_erf_inputs.py --ranks 64 --umax 30 Exec/RegTests/*/inputs*
python3 Exec/Tools/check_erf_inputs.py --min-severity WARNING --strict "Tests/test_files/*/*.i"
python3 Exec/Tools/check_erf_inputs.py --json inputs > report.json
python3 Exec/Tools/check_erf_inputs.py --self-test
```

| Option | Default | Effect |
|--------|---------|--------|
| `--ranks N` | unset | box-count balance advice for N MPI ranks |
| `--gpu` | off | box-size advice for GPUs (larger boxes) |
| `--umax U` | 20 m/s | wind speed used for the CFL estimates |
| `--theta T` | 300 K | potential temperature for the sound speed |
| `--min-severity L` | INFO | hide findings below L |
| `--strict` | off | exit 1 on WARNING as well as ERROR |
| `--only NAME` | all | run one checker (`--list-checkers` lists them) |
| `--include-bases` | off | also check files that a sibling deck includes with `FILE =`, and `*master*` decks |
| `--docs PATH` | auto | `Inputs.rst` location when run outside the repository |

Run it from inside the repository so the Sphinx pages and `Source/` are
found; the misspelt-key check is skipped otherwise.

## Real-data and external files

When the deck names a file that is present, its contents are compared with the
deck. wrfinput: staggered dimensions against `amr.n_cell`, `DX * n` against the
domain extent, `bottom_top` against `n_cell` z, `SIMULATION_START_DATE` against
`start_datetime`. wrfbdy: `bdy_width` above `erf.real_width`, evenly spaced
time stamps, data span against `stop_time`, the cloud fields that
`erf.use_wrf_bdy_qc_qi` needs. wrflow: SST present, grid against `n_cell`.
Terrain text, USGS and NetCDF files: shape, increasing coordinates, coverage of
the domain (ERF stops when the table ends more than one spacing inside the
domain). Soundings: columns, increasing heights, top against the domain top,
surface theta against `erf.most.surf_temp`. NetCDF files are read with the
`netCDF4` module when installed, otherwise with `ncdump -h`; without either
only the existence of the file is checked.

## Expected results

A clean deck prints INFO lines only and the script exits 0. The self-test
exits 0 with `self-test: 109 cases passed`; CTest runs it as
`InputsChecker_SelfTest`. Over the repository's decks (`Exec/RegTests`,
`Tests/test_files`, `Exec/CanonicalTests`, 729 files) the first run on
2026-10-23 reported 24 decks with errors, most of them external files that the
test setup provides later; the remaining groups are filed as
erf-model/ERF#4424 to #4429.

## Adding a check

Add a class in `checks/`, decorate it with `@register`, set `name`, implement
`run()` and, for a mode-specific class, `applies(ctx)`. Use `self.error`,
`self.warn`, `self.advise`, `self.info` with the key, what is wrong and what
to do. Add one tripping and one clean case to `selftest.py`; the self-test
fails when either does not behave.
