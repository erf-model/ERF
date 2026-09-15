# ERF SBM P2 source trace

This is the implementation trace for the P2 milestone.  The historical
P0/P1 trace remains in `P0_P1_SOURCE_TRACE.md`; this file records the P2
contracts and their actual ERF/AMReX integration points.

## G0 — source contracts

| Contract | Actual source symbol | Source path | Units / semantics | P2 use | Unsupported cases |
|---|---|---|---|---|---|
| Carrier mass flux | `ERF::advance_sbm_stage`, `avg_xmom/avg_ymom/avg_zmom` | `Source/Microphysics/SBM/ERF_SBMErfIntegration.cpp`, `Source/ERF.H` | ERF face-centered dry-air mass flux; the transport kernel consumes the existing face arrays | Common carrier flux for donor and WENO candidates | Independently reconstructed velocity fluxes |
| High-order reconstruction | `WENO_Z3::InterpolateInX/Y/Z` | `Source/Utils/ERF_Interpolation_WENO_Z.H`; ordinary scalar caller `Source/Advection/ERF_AdvectionSrcForState.cpp` | Reusable ERF scalar WENO-Z3 helper on a grown `Array4`; reconstructed quantity is the intensive `X/rho` | `GroupedFCT_WENOZ3` candidate | Non-orthogonal/terrain reconstruction |
| Auxiliary stage state | `AuxiliaryStateManager::{old,evaluation,output,scratch}` | `Source/AuxiliaryState/ERF_AuxiliaryStateManager.{H,cpp}` | Old full-step baseline, accepted evaluation/output, bounded scratch | Stage recurrence and 2M endpoint workspace | Permanent full endpoint state |
| Stage timing | `make_compressible_stage`, `make_anelastic_stage`, `StageContext::accepted_ledger_weight` | `Source/AuxiliaryState/ERF_AuxiliaryStageContext.H` | Compressible callback is anchored at `Z^n`; Heun corrector uses half-step transfer | Low/high update, FCT correction, accepted ledger | Final-stage-only Heun ledger |
| Physical face transfer | `AuxiliaryFaceTransferLedger::record_stage` | `Source/AuxiliaryState/ERF_AuxiliaryFaceTransfer.{H,cpp}` | Conceptual `I=A*integral(F dt)`; internal face FAB is per-area flux and stage weights are explicit | One accepted transfer feeds divergence, projection, diagnostics, and AMR | Candidate/pre-limit transfers |
| Diffusivity | P2 manufactured coefficient `solverChoice.sbm_diffusion_coeff` | `Source/DataStructs/ERF_DataStruct.H`, `Source/Microphysics/SBM/ERF_SBMDiffusion.cpp` | Explicit coefficient `K`; face operator is `-rho_f*K*grad(X/rho)` | Orthogonal two-point diffusion in the same accepted flux | Native implicit moisture diffusion, SHOC, EB, terrain cross terms |
| Geometry | `Geometry::CellSizeArray`, `Geometry::InvCellSize`, `Geometry::periodicity` | `Source/Microphysics/SBM/ERF_SBMTransportPrototype.cpp` and ERF level hooks | Cartesian face flux with divergence factors `1/dx`, `1/dy`, `1/dz`; generic adapter tests use area and volume explicitly | Physical transfer scaling and timestep bound | Moving terrain and EB metrics |
| Flux register scaling | `YAFluxRegisterT::CrseAdd/FineAdd` | pinned `Submodules/amrex/Src/Boundary/AMReX_YAFluxRegister.H`; ERF adapter in `ERF_SBMErfIntegration.cpp` | API expects per-area instantaneous flux and multiplies by supplied `dt/dx`; P2 passes accepted stage flux and weighted `dt` exactly once | Provider-owned `sbm_flux_reg` | Area-integrated data passed without conversion |
| FillPatch/time interpolation | `ERF::FillCoarsePatch`, `FillPatchFineLevel` | `Source/BoundaryConditions/ERF_FillCoarsePatch.cpp`, `Source/BoundaryConditions/ERF_FillPatch.cpp` | Core ERF hierarchy path; P2 auxiliary creation uses conservative coarse injection and periodic spectral ghost fill | Fine-level auxiliary lifecycle | Extrapolated auxiliary ghosts |
| Restriction/average-down | `ERF::AverageDown`, `AverageDownTo`, `AuxiliaryStateManager::average_down_to` | `Source/Utils/ERF_AverageDown.cpp`, `Source/AuxiliaryState/ERF_AuxiliaryStateManager.cpp` | AMReX volume-consistent average-down of authoritative components | Compact fields are projected after spectral average-down | Independent `qc/qr` averaging |
| Level creation | `ERF::MakeNewLevelFromScratch`, `MakeNewLevelFromCoarse` | `Source/ERF_MakeNewLevel.cpp` | Auxiliary allocation at every SBM level; coarse-to-fine uses `pc_interp` and then projects compact state | Fresh and refined levels | Compact-to-spectrum reconstruction |
| Regrid lifecycle | `ERF::RemakeLevel`, `ERF::ClearLevel` | `Source/ERF_MakeNewLevel.cpp` | Remake copies state into new FABs and rebuilds face ledgers; clear destroys all per-level ownership | No stale FAB/device views after regrid | Stale register/pointer retention |
| Checkpoint/restart | `ERF::WriteCheckpointFile`, `ERF::ReadCheckpointFile` | `Source/IO/ERF_Checkpoint.cpp`; schema service `ERF_SBMRestart.{H,cpp}` | Per-level `SBMAux`, exact schema, compact projection compared before overwrite | Strict P2 restart integrity | Schema conversion or bulk reconstruction |
| Physical boundaries | `BoundaryKind`, ERF `phys_bc_type` classification | `Source/Microphysics/SBM/ERF_SBMBoundary.{H,cpp}`, `Source/BoundaryConditions` | Periodic, prescribed spectral inflow, advective outflow, impermeable wall; budgets are signed `A*dt*F` | Generic boundary transfer service and tests | Wall deposition, unspecified incoming outflow spectrum |

## G1 — layout and invariant domain

`SBMLayout` resolves runtime population/bin/property offsets once on the host.
Checkpoint/storage is one-moment `M` or two-moment physical `(M,C)`; the
two-moment path uses the same bounded scratch slots for endpoint `(L,H)` only
during transport.  `make_constraint_groups` emits one complete group per
population/bin and linear forms for mass, endpoints, attached properties,
support bounds, and mass-bounded subsets.  `transform_two_moment` uses FMA and
the cancellation-scaled tolerance `128*epsilon*scale`, normalizing only tiny
negative endpoints.

The P2 tests cover 4/16/64 runtime bins, exact cone edges, empty states,
support/subset constraints, physical 2M storage, and a two-population
synthetic ice-like subset property.  No ice process is present.

## G2 — accepted grouped transport

`ERF_SBMTransportPrototype::advance_stage` constructs donor fluxes once from
the existing ERF carrier mass flux.  The grouped path fills bounded scratch
with `X/rho`, calls the shared ERF `WENO_Z3` helper, and applies one face
limiter across the complete bin group.  The host `limit_grouped` implementation
is the deterministic reference for global cell/face/constraint budgets; the
production FAB path uses POD property descriptors and the same one-lambda
group semantics.  Final accepted physical `(M,C)` transfers are recorded and
projected; no compact field has an independent transport path.

The P2 test executable exercises positive/negative transfer signs, common
limiting, duplicate face-ownership rejection, chunk-policy invariance, the
production WENO path, and both temporal stage contracts.

## G3 — diffusion and boundaries

`ERF_SBMDiffusion` documents and tests the physical two-point operator, while
the production transport path adds the same density-weighted orthogonal term
to its face flux.  The production timestep is rejected when the explicit
diffusion bound is exceeded; there is no hidden auxiliary subcycling or
clipping.  The boundary service rejects incomplete/nonrealizable prescribed
inflow and emits zero transfer for periodic and impermeable-wall descriptors.

## G4 — AMR and accounting

`sbm_flux_reg` is allocated only when SBM and two-way coupling are active.  The
integration hook registers the accepted spectral face FAB using YAFluxRegister's
per-area-flux plus `dt/dx` convention.  `post_timestep` refluxes the
authoritative spectrum, validates it, projects `qc/qr`, then performs the
normal ERF average-down and matched auxiliary average-down.  An inadmissible
post-reflux state is reported by `validate_post_reflux` and fails closed; no
repair/clipping path exists.  The manager lifecycle test covers create,
prolong, average-down, remake with a changed box decomposition, and destroy.
The `SBM_P2_AMR_2M` production fixture additionally runs the grouped-FCT 2M
path on two levels and two MPI ranks through coarse step 2.

## G5 — restart and regrid

`SBM_Schema` records layout, moment, grid, property, projection, transport,
constraint, and numerical identities.  On restart, ERF reads and compares the
schema and per-level authoritative `SBMAux`, reconstructs the compact
projection in scratch, compares it against checkpointed compact `qc/qr`, and
only then overwrites the normal core fields.  Changed schema fields are not
converted automatically.

## G6 — qualification boundary

The P2 capability gate is selected by `sbm_transport_method`, 2M mode, or a
positive explicit SBM diffusion coefficient.  It fails closed for SHOC,
implicit moisture diffusion, moving terrain, EB, nonperiodic production
geometry, dynamic grids, schema conversion, and P3+ physical processes.
Ordinary ERF paths remain unmodified when SBM is disabled.  Qualification
counts and environment limitations are recorded in `P2_QUALIFICATION_REPORT.md`.
