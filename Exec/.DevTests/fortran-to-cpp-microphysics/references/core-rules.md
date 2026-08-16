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
- **The source is the oracle for STRUCTURE AND PLACEMENT, not only for values.**
  Before deciding where a computation lives, what gates it, or how a loop is
  split, locate the corresponding source construct and quote it. Answer three
  questions explicitly:
  - *Where* does the source do this — which routine, and at what point in the
    per-step sequence?
  - *What* does the source gate it on — the actual condition, not a proxy that
    happens to agree in the common case?
  - *How* is it nested — one loop or several, and which values are computed once
    per cell and reused across sub-blocks?

  This is cheap, usually one grep, and it is the highest-value check available
  because a structural error is invisible to value comparison until it fires.
  Four distinct defects in the WDM6 campaign traced to skipping it: seeding
  placed in `Init` when the source does it in the per-step driver; a restart
  string used as a proxy for the source's step condition; a per-cell scratch
  value proposed where the source shares one loop nest; and a sibling scheme
  assumed to have a precedent it did not have. Each was confidently wrong rather
  than uncertain, so more deliberation would not have caught any of them — only
  the check does.

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
- Route **every** unsuffixed source literal through the single conversion
  helper, including inline literals inside kernels, not just named parameters.
  Filter attention to literals that are not binary-exact and are used in
  arithmetic; integers and dyadic fractions such as `2.`, `4.`, `24.`, `0.5`
  and `1.e3` need no routing. Prioritize **exponents**, which amplify: the
  error enters as `ln(base) * delta`, so a 4e-08 exponent difference became
  4.3e-07 on the result, and then 2.1e-06 after a downstream subtractive
  cancellation.
- **A literal inside a clamp or a gate is a NULL RESULT until something
  reaches it, so an early clean milestone says nothing about it.** The WDM6
  rain-slope cap, `min(1./lamdar(...),1.e-3)` at `ERF_module_mp_wdm6.F90:3420`
  and `:3493`, was unrouted on the C++ side for the whole campaign and cost
  nothing until step 16, because before that no cell's slope reached the cap.
  Milestone B was bitwise-clean through all ten of its steps *and stayed
  numerically identical after the fix* — the invariance is what proves the
  clamp never bound there. Do not read an early milestone pass, or a fix that
  leaves it unmoved, as evidence about a clamped or gated literal. Audit these
  by reading the source, not by watching a metric.
- The audit is not only a defect hunt: **it is the precondition for the
  double-literal comparison build to mean anything.** An unrouted literal is
  already a true double on the C++ side, so it accidentally agrees with the
  Fortran under the double-literal configuration and disagrees only under the
  float configuration. Every unrouted literal therefore silently removes the
  discrimination that build exists to provide, and a clean comparison across
  the two modes can indicate uniform routing rather than genuine agreement.
  Complete the routing audit before drawing any conclusion from a mode-to-mode
  comparison.
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
