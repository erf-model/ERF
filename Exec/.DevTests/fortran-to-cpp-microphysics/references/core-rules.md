# Core Rules

These are the reusable, cross-scheme rules that should survive beyond WSM6.

## Vocabulary
- `path_a_fortran_bridge`: tracked Fortran implementation wired into ERF and used as the validation oracle.
- `path_b_native_cpp`: native AMReX C++ implementation validated against Path A.
- `phase_1_bridge` through `phase_5_full_native_validation`: the canonical staged workflow from bridge wiring to plotfile parity.
- `tag_frontier`: controlled-scope retreat lane.
- `plotfile_parity`: milestone-based full-domain parity lane.

## Bridge and Interface Rules
- Register new schemes through the full ERF moisture-model path: `MoistureType`, `modelType()`, `EulerianMicrophysics`, and `Source/Microphysics/<Scheme>/Make.package`.
- Keep the original Fortran source unchanged. Put adaptation logic in `iso_c_binding` wrappers and C interfaces.
- Gate bridge-vs-native behavior through runtime toggles when the scheme already supports both paths. Do not force users to rebuild just to compare implementations.
- Pass storage bounds and active bounds separately when the Fortran scheme distinguishes memory extents from compute extents.

## Native Port Rules
- Use the two-layer state model:
  - persistent `MicVar_<Scheme>` storage owned by ERF
  - working `<Scheme>Ind` storage that mirrors the Fortran kernel
- Keep pack/unpack boundaries explicit so Path A and Path B remain auditable.
- Resolve physical constants via an explicit hierarchy:
  - ERF constants first
  - class members for derived or variant-controlled values
  - in-kernel computation for temperature- or state-dependent quantities
- Keep init-time coefficient derivation logically separate from runtime advance kernels.
- Choose one of the established AMReX translation strategies and stay consistent:
  - enum + working FAB
  - device struct with inline methods
- Survey naming conventions before minting new enum entries or persistent state names.
- Translate intent, not dead Fortran scaffolding. Omit dead code paths and inline truly trivial helpers only when parity is preserved.

## Numerical Fidelity Rules
- When the source Fortran writes unsuffixed literals, the float assumption
  applies to the whole literal **expression**, not to each literal separately.
  `.104*9.8` is a default-REAL expression: the product is evaluated in single
  precision and only then widened. A port that converts each literal
  individually and multiplies in double does not reproduce it. Convert the
  operands, then apply the same conversion to their product; do not instead
  convert the double-precision product, which is a different and generally
  worse value. Check for shadowed constants on both sides of the port while
  reading these declarations.
- Port init-time special functions exactly when they influence later parity. Approximate replacements are not acceptable unless separately validated.
- Treat helper semantics as part of parity, not as harmless cleanup.
- Respect slab/box allocation patterns for 2D accumulators and other mixed-dimensional storage.
- Perform a pre-compile alignment pass group-by-group to catch dropped boundary statements and loop-bound mistakes before debugging runtime behavior.
- Explicitly restructure risky Fortran control-flow idioms instead of transliterating them blindly.

## Large-Routine Workflow Rules
- Use the staged five-phase workflow for large ports.
- Build a complete process inventory before writing the native kernel.
- Treat minor-timestep outer loops and column-serial sedimentation kernels as generic bulk-microphysics concepts, even when WSM6 is the demonstrated example.

## Extraction Policy
- Promote only rules that generalize cleanly beyond a single scheme or test case.
- Keep scheme-local execution orderings, tag tables, milestone step numbers, and variable-to-tag mappings out of this file unless multiple schemes confirm them.
