# erf_inputs_schema

The package behind `Exec/Tools/build_erf_schema.py`: a per-build schema of ERF's inputs, read from the source of one
ERF commit and its AMReX. It runs offline and read-only. Source is read by `git archive`, so nothing is checked out
and nothing is built.

| module | does |
|---|---|
| `scan.py` | finds every literal ParmParse query (`query`, `get`, `queryarr`, `getarr`, the enum forms, `query_one_or_per_level*`, `read_int_string`, `contains`/`countval`), with the C++ type of the variable it fills and the `#ifdef ERF_USE_X` it sits under |
| `cmake.py` | which `ERF_USE_*` macros a build defines: ERF's option defaults, the build's `-D` options or `CMakeCache.txt` over them, and the `if()` around each `target_compile_definitions` |
| `schema.py` | merges the reads into one fact per key and assembles the schema |
| `canonical.py` | the type codes and the canonical form of each value |
| `source.py` | `git archive` of an ERF commit and of AMReX at ERF's `Submodules/AMReX` gitlink |

## One key's fact
| field | meaning |
|---|---|
| `type` | `int`, `real`, `bool`, `string`, `enum`, each also as `_list`; `unknown` (the variable's type was not found); or `removed` |
| `source`, `layer`, `layers` | the read the fact is taken from; `AMReX` or `ERF`; and each layer's first read (a key AMReX and ERF both read keeps both lines) |
| `per_level` | one value, or one per level |
| `values`, `case_insensitive` | an enum's names, and whether ERF reads them case-insensitively |
| `default`, `default_origin`, `default_source` | the variable's initializer as a canonical literal (`initializer`), `none`, or `expression` when it is not a plain literal |
| `required` | some read is a `get` |
| `requires`, `requires_evidence` | the `ERF_USE_X` macros every read sits under, with the `#ifdef` line |
| `value_requires` | an enum value (or a site fact's value) that needs a feature, with evidence |
| `removed_shape`, `removed_evidence`, `shadowed_reads` | for `removed`: the shape (`contains_abort`, `query_gt0_error`, `countval_gt0_error`, `countval_var_gt0_error`), the Abort/Error line, and any typed read elsewhere that the removal overrides |
| `type_conflict`, `family` | reads that disagree on the type; a family pattern (`<face>`, `<station>`, `<n>`) |

**A removed key.** A key is `removed` when its read sits in an `if()` whose body always calls `Abort` or `Error`, at
the top level of the body, and the body does not go on to read the key. A size check nested in the body
(`terrain_z_levels`) is not a removal, and neither is a deprecated key that is read and mapped (AMReX's `nosub`).

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
