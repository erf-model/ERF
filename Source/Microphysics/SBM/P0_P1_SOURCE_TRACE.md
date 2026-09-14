# ERF SBM P0/P1 source trace

This trace records the source locations and contracts used by the experimental
P0/P1 implementation.  It is an implementation map, not a replacement for
the design document or a claim of physical microphysics qualification.

## Provenance

* P1 base: `origin/sbm-p0-p1-aux-state` at `c321be7007da476e4218186fa76698f7f22e02e8`
* P1 branch: `sbm-p1-qualification`
* ERF development ancestor: `33ce039e87f309592609098792dff06d0762438c`
* AMReX submodule: `e60cdc18711ccf7fc0616d7a2fdf062021376976`
* Qualification prompt SHA-256: `ffd655e0087e72c5d1b125100d46fabf8d7b88e19848bb8aadcbb29e9d5c688`
* Final hardening prompt SHA-256: `ae5adcca287e8df24de2c2910c50ccad4318bc988b8ecee19691904b54d9451e`
* Final qualification fixes prompt SHA-256: `cfbfa4ddbdbcdff30ad05b4d92190923b52dc137fd3c862bbd157586ff907f17`

## Authoritative state and stage path

`Source/AuxiliaryState/ERF_AuxiliaryStateManager.{H,cpp}` owns generic old,
evaluation, output, scratch, and face-ledger storage.  The stage algebra is
explicit in `ERF_AuxiliaryStageContext.H`: compressible ERF RK3 uses the old
full-step baseline at every callback, while anelastic Heun uses the predictor
only as the stage-1 evaluation and corrects against the unchanged old state.

`Source/TimeIntegration/ERF_TI_slow_rhs_post.H` calls
`ERF::advance_sbm_stage` after ERF's host slow-RHS, state update, and boundary
operations.  `Source/Microphysics/SBM/ERF_SBMErfIntegration.cpp` constructs the
provider context and passes ERF's actual `avg_xmom`, `avg_ymom`, and `avg_zmom`
carrier fields to `ERF_SBMTransportPrototype.cpp`.

`ERF::Advance` swaps ERF's old/new state vectors before entering the dycore.
`ERF::begin_sbm_step` therefore receives the explicit `state_old[cons]`
MultiFab from `advance_dycore` and snapshots its compact `qc`/`qr` components.
This keeps the compact baseline tied to the true start of each full step,
including the second and subsequent steps; it does not read the destination
state after the swap.

## Actual production face-transfer ledger

`Source/AuxiliaryState/ERF_AuxiliaryFaceTransfer.{H,cpp}` provides compact,
generic stage-local x/y/z face-centered `MultiFab` storage and
`AuxiliaryFaceTransferLedger`.  The ledger owns one reusable stage scratch
transfer and one accepted accumulator; it does not retain a three-stage
history.  Each call to `record_stage` copies the exact stage flux into the
scratch object and accumulates the accepted transfer using `StageContext`:

* compressible: `I = dt * F(stage 2)`;
* anelastic: `I = dt/2 * (F(stage 0) + F(stage 1))`.

`Source/Microphysics/SBM/ERF_SBMTransportPrototype.cpp` constructs each donor
face flux once in the stage transfer, uses those same arrays in the spectral
divergence, validates the updated authoritative spectrum, and records the
stage.  `ERF_SBMBulkProjection.cpp` projects the accepted spectral face
transfer into two compact bulk face components.  No independently reconstructed
`qc`/`qr` face flux is used.  The accepted face object is retained in the ERF
provider as `sbm_accepted_bulk_face_transfer`, ready for a future AMR reflux
consumer; P1 itself remains single-level.

`Source/Microphysics/SBM/ERF_SBMTransferClosure.{H,cpp}` evaluates the local
accepted-transfer identity `X_new - X_old + div(I_accepted)` for every
spectral component and for compact `qc`/`qr`.  The diagnostic stores the
maximum residual and its machine-epsilon-scaled tolerance.  The standalone
`Tests/SBMQualificationCheck.cpp` recomputes the pass conditions from those
values independently of the producer's `passed` flag.

## Compact projection and metadata

`ERF_SBMLayout.{H,cpp}` separates generic `SpectralGrid` coordinate metadata
from liquid-only projection metadata.  A population has a generic
`mass_offset`; `LiquidProjectionSpec` owns the cloud/rain split.  A second
non-liquid population can therefore be represented without liquid semantics.
`CoordinateKind::{Mass,Radius}` and the corresponding units are phase-neutral;
the schema identities are `spectral-grid-v3` and `sbm-layout-v2`.  Attached
properties explicitly declare `ExtensiveMass`, `NumberCarried`, or
`MassBoundedSubset` semantics and their carrier-bin-conservative remap policy.
The grid coordinate is in `kg` for the qualified liquid/aerosol cases, while
transported mass components are in `kg m^-3` and number components are in
`m^-3`.

`SBMBulkProjection::apply_to_core` is the only P1 compact-state projection:
`qc` is the sum of the liquid cloud bins and `qr` is the sum of the liquid rain
bins.  `apply_to_face_transfer` applies the identical partition to accepted
face transfers.

## Ownership and capability gating

`ERF_SBMContracts.{H,cpp}` contains both the stable capability report and the
host-side `OwnershipRegistry`.  `ERF_SlowRhsPost.cpp` asks that registry for
the number of native moisture components, leaving `qv` on ERF's path while
`qc`/`qr` are provider-owned.  `ERF_DataStruct.H` invokes the same capability
evaluator for fail-closed configuration validation; unsupported diffusion,
forcing, terrain, coupling, P2 physics, and two-moment transport are rejected
before the run proceeds.  `ERF_SBMErfIntegration.cpp` separately rejects
restart/schema conversion before SBM state initialization.  No strings or
dynamic ownership lookups enter device kernels.

## Runtime inputs and finite checks

`ERF_DataStruct.H` reads `sbm_*` settings only after selecting the SBM
moisture model.  `sbm_nbins` and the cloud/rain split are validated before
allocation; edge/pivot lengths, finiteness, monotonicity, and bin membership
are checked.  `validate_runtime_bin_count` is the small pure contract used by
the parser and unit tests.

`ERF_SBMTransportPrototype::validate_nonnegative_state` performs one fused,
GPU-capable `ParReduce` over all requested state components, returning a
nonfinite flag, minimum, and maximum absolute value.  Because AMReX's
`ParReduce` is local, exactly three scalar MPI reductions combine those values
globally, independent of the number of bins.  NaN, positive/negative infinity,
and material negative values fail closed.  The negative tolerance is
`128 * epsilon * max_abs(state)` with no order-one floor, so an all-zero state
has zero tolerance; values within that scale are accepted without clipping.

The diagnostic memory model reports cell-state bytes, face-transfer bytes, and
their total as logical allocated data payloads.  It counts each grown FAB box,
including ghost cells, for the four cell states, the reusable stage plus
accepted spectral ledger objects, the accepted compact bulk projection, and
the initial compact state snapshot.  Allocator metadata is excluded.  The
versioned diagnostic also records the completed full-step count and the real
qualification cases require at least two.

## P1 capability boundary

Qualified infrastructure: one-moment liquid mass transport, runtime bin
counts, static single-level periodic Cartesian manufactured transport,
compressible RK3, anelastic Heun, exact accepted face-transfer projection, and
compact projection conservation.

Explicitly unsupported: AMR execution/reflux, moving or non-Cartesian terrain,
embedded boundaries, diffusion, implicit moisture diffusion, SHOC/macrophysics,
high-order/FCT transport, sedimentation, condensation/evaporation,
activation/regeneration, collision/coalescence, ice, aerosol lifecycle physics,
dynamic grids, two-moment transport, independent moisture forcing,
large-scale/sounding forcing, sponge/wall modification, restart/schema
conversion, GPU performance qualification, and physical warm-cloud validation.

## Tests

* `Tests/Unit/Microphysics/SBM/ERF_GTestSBMP0P1.cpp`: grid/layout metadata,
  two-moment algebra, temporal recurrence, face-ledger weights, topology,
  ownership, finite/negative controls, accepted-state local closure (including
  wrong-final-stage, perturbed-face, zero-scale, and stale-baseline controls),
  free-stream preservation, and runtime 4/16/64-bin transport.
* `Tests/CTestList.cmake` and `Tests/RunSBMPrototype.cmake`: six real ERF MPI
  cases (compressible/anelastic × 4/16/64 bins).
* `Tests/SBMQualificationCheck.cpp`: independent numerical diagnostic checker;
  it verifies the emitted invariant values so CTest does not rely on log text
  alone.

The six manufactured cases use a nonuniform periodic spectrum, run for two
complete level-0 steps, verify spectral mass conservation and compact
projection, and check the accepted face ledger.  They emit the versioned
diagnostic's step count, three local closure residuals/tolerances, and the
three-part grown-FAB memory accounting consumed by the independent checker.
The documentation
build was attempted with `cd Docs && ./BuildDocs.sh`; its Python catalog checks
pass, but this checkout has neither `doxygen` nor `sphinx-build` (and
`python3 -m sphinx` is unavailable).  Reproduce that exact attempt after
installing those documented tools.  The existing variable-density free-stream
test remains a separate preservation test rather than a substitute for
nonuniform transport.
