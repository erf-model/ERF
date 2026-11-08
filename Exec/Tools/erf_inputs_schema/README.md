# erf_inputs_schema

The package behind `Exec/Tools/build_erf_schema.py`: a per-build schema of ERF's inputs, read from the source of one
ERF commit and its AMReX. It runs offline and read-only. Source is read by `git archive`, so nothing is checked out
and nothing is built.

| module | does |
|---|---|
| `scan.py` | finds every literal ParmParse query (every read method of `AMReX_ParmParse.H`: `query`/`get` and their `arr`, `kth`, `ktharr`, `line`, `table`, `WithParser`, `AsDouble`, `Add` and enum (`_case_insensitive`, `_sloppy`) forms; `query_one_or_per_level*`, `read_int_string`, `contains`/`countval`), with the C++ type of the variable it fills and the `#ifdef ERF_USE_X` it sits under; and the literal `pp.add` calls of ERF's `add_par` |
| `cmake.py` | which `ERF_USE_*` macros a build defines: ERF's option defaults, the build's `-D` options or `CMakeCache.txt` over them, and the `if()` around each `target_compile_definitions` |
| `schema.py` | merges the reads into one fact per key and assembles the schema |
| `canonical.py` | the type codes and the canonical form of each value |
| `source.py` | `git archive` of an ERF commit and of AMReX at ERF's `Submodules/AMReX` gitlink |

## One key's fact
| field | meaning |
|---|---|
| `type` | `int`, `real`, `bool`, `string`, `enum`, each also as `_list`; `unknown` (the variable's type was not found), `unknown_list` (an array read of such a variable); or `removed` |
| `source`, `layer`, `layers` | the read the fact is taken from; `AMReX` or `ERF`; and each layer's first read (a key AMReX and ERF both read keeps both lines) |
| `per_level` | one value, or one per level |
| `values`, `case_insensitive`, `ignores` | an enum's names, whether ERF reads them case-insensitively, and the characters a `query_enum_sloppy` read drops from the word first |
| `line`, `table` | read with `queryline` (the whole line, any number of words) or `querytable` (a `{{a, b} {c, d}}` table, typed `unknown`) |
| `parser`, `as_double` | read with a `WithParser` or `AsDouble` method, which takes an expression by design; `AsDouble` also reads an int as a double and rounds it |
| `default`, `default_origin`, `default_source` | the variable's initializer as a canonical literal (`initializer`), `none`, or `expression` when it is not a plain literal; `add` when ERF's `add_par` sets it with `pp.add("k", literal)`, which runs before the inputs are read, so the inputs override it (`amr.n_proper` is 2, not AMReX's 1) |
| `required` | some read is a `get` |
| `requires`, `requires_evidence` | the `ERF_USE_X` macros every read sits under, with the `#ifdef` line |
| `value_requires` | an enum value (or a site fact's value) that needs a feature, with evidence |
| `removed_shape`, `removed_evidence`, `shadowed_reads` | for `removed`: the shape (`contains_abort`, `query_gt0_error`, `countval_gt0_error`, `countval_var_gt0_error`), the Abort/Error line, and any typed read elsewhere that the removal overrides |
| `type_conflict`, `family` | reads that disagree on the type; a family pattern (`<face>`, `<station>`, `<indicator>`, `<n>`) |
| `refused_by` | a refusal of the key in problem-specific code (`Source/Prob/`), which does not override the typed read elsewhere |

**A removed key.** A key is `removed` when its read sits in an `if()` whose body always calls `Abort` or `Error`, at
the top level of the body, and the body does not go on to read the key. A size check nested in the body
(`terrain_z_levels`) is not a removal, and neither is a deprecated key that is read and mapped (AMReX's `nosub`).

**Two reads that disagree.** `unknown` yields to any type; a `queryarr`/`getarr` read of a key makes it a list even when that read's variable is untyped (`erf.vert_implicit_fac`); a scalar read of one element (`fixed_dt[0]`) is a scalar; an enum read wins over a `std::string` probe of the same key (`erf.init_type`, `erf.moisture_model`, `erf.terrain_type`), with the pair kept in `type_conflict`. A default is taken only from a read of the type the key ends with, and never from text inside a string literal.

**Prefixes built at run time.** `scan.FAMILIES` maps a `ParmParse` constructor argument in one file to its prefixes: the boundary faces (`<face>`), stations (`<station>`), refinement indicators (`erf.<indicator>`), obs-nudging stations, the particle containers (`tracer_particles`, `super_droplets_moisture`). A `ParmParse&` parameter takes the prefixes of the objects passed at its call sites (`parse_wall_face_scalars`). The generator stops when an entry's file exists but nothing matches it, so an edited constructor argument is caught. Reads that still do not resolve are listed under `unresolved`, and the checker names them instead of calling the key unread.

**One name with several kinds in a file.** The kinds the read call takes are kept: a list for the array reads, a
scalar otherwise, any for a per-level read. If several remain, the nearest declaration before the read decides.
Otherwise the type is `unknown`.

## Build facts and site facts
- **Build options:** `--cmake='-DERF_ENABLE_X=ON …'` (repeatable; use the `=` form, because the value starts with a
  dash) and/or `--cmake-cache CMakeCache.txt`. A cache holds every ERF option. A `-D` list holds only the options
  passed, and the rest take the CMakeLists defaults. `build.facts_from.cmake` records which kind was used.
- **Site facts:** `--site-facts FILE` is a JSON list of facts the build options cannot show, for example:
  ```json
  [{"feature": "mpi_gpu_aware", "on": true, "key": "amrex.use_gpu_aware_mpi", "value": "true",
    "evidence": "run time: MPICH_GPU_SUPPORT_ENABLED=1"}]
  ```
  Each one sets `features[feature]`. When it names a key and a value, it also sets that value's requirement.

## Metadata and determinism
- **Recorded:** `schema_version`, `erf_commit`, `amrex_commit`, `generator`, `generator_commit`, `build`
  (`cmake`, `facts_from`) and `canonical_rules`.
- **Determinism:** the same commits, options, site facts and generator commit give the same bytes (sorted keys, no
  timestamps).

## Tests
```
python3 -m unittest discover -s Exec/Tools/erf_inputs_schema/tests
python3 Exec/Tools/check_erf_inputs.py --self-test
```
The fixture is a small synthetic source tree (`tests/fixtures/synthetic/`) whose files copy ERF's own source shapes,
each cited in a comment. `synthetic_expected_keys.json` is its schema, written by hand.
