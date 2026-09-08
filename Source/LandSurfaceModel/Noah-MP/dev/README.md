# ERF Noah-MP driver — developer specs (`dev/`)

This directory holds the design specifications for the **ERF-side driver of
Noah-MP**: the `NOAHMP` C++ class that plugs the Noah-MP land-surface model into
ERF's land-surface-model interface and drives the per-step state exchange between
the two. The driver is split by concern across several files:

| File                     | Role                                                         |
|--------------------------|--------------------------------------------------------------|
| `ERF_NOAHMP.H`           | The `NOAHMP` class declaration: members, accessors, helpers. |
| `ERF_NOAHMP_Fields.H`    | The coupling-field enums and the X-macro field registry.     |
| `ERF_NOAHMP_Init.cpp`    | `Init` and the run lifecycle.                                |
| `ERF_NOAHMP_Advance.cpp` | `Advance_With_State` and the per-step pipeline helpers.      |
| `ERF_NOAHMP_Precip.cpp`  | Precip source collection, snapshots, and guard reporting.    |
| `ERF_NOAHMP_IO.cpp`      | The land plotfile and checkpoint/restart.                    |

These are *living design documents* for contributors, not a user manual. They
describe what the driver does, why it is built the way it is, and the invariants
a change must not break. Public, user-facing documentation lives at
<https://erf.readthedocs.io/en/latest/CouplingToNoahMP.html>.

> **Scope boundary.** This driver is the *ERF half* of the coupling. The layer
> that lets C++ call into the Fortran land model (the auto-generated language
> bindings) and the Noah-MP physics itself live in the Noah-MP submodule under
> `Submodules/Noah-MP/drivers/erf/` and have their own `dev/` specs. When a topic
> crosses into that layer (the `NoahmpIO_type` API, array index ordering, NetCDF
> restart format), these specs link there rather than restate it.

## Index

| Spec                                                                  | What it covers                                                                                                                                                                  |
|-----------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [<code>spec&#8209;noahmp&#8209;api.md</code>](spec-noahmp-api.md)     | The `NOAHMP` class, the `NullSurf` contract, the `LsmData`/`LsmFlux` coupling fields, the run lifecycle and time subcycling, and adding a coupled variable. **Start here.**     |
| [<code>spec&#8209;noahmp&#8209;gpu.md</code>](spec-noahmp-gpu.md)     | The GPU-aware, component-indexed state exchange: staging buffers, the `ParallelFor`/`LoopOnCpu` dataflow, sync points, slab/ghost-cell handling, and the flux fill-value guard. |
| [<code>spec&#8209;noahmp&#8209;io.md</code>](spec-noahmp-io.md)       | Driver-owned I/O: checkpoint/restart, the per-step land plotfile, and static-input reads (`namelist.erf`, `NoahmpTable.TBL`, the NetCDF land file).                             |
| [<code>spec&#8209;noahmp&#8209;reorg.md</code>](spec-noahmp-reorg.md) | The driver's **source layout** and the X-macro field registry: why the code is split as it is, and how the enum/name/copy sync is collapsed.                                    |

## Conventions used across these docs

- **ERF side vs Noah-MP side.** "ERF side" = ERF's gridded fields (held in AMReX
  `MultiFab`/`FArrayBox` containers), indexed by grid cell `(i,j,k)` with the
  ground at the lowest vertical level `k = klo`. "Noah-MP side" = the
  `NoahmpIO_type` Fortran-backed arrays, indexed in Fortran's order with the
  **vertical `k` and second horizontal `j` axes swapped** — atmospheric fields are
  `(i,k,j)` with `k` pinned to the surface layer `1`. The driver's job is to move
  data between these two layouts correctly; see [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md).
- **Coupling fields.** The fields exchanged with the rest of ERF (radiation,
  surface layer) are listed in `LsmData_NOAHMP` (state) and `LsmFlux_NOAHMP`
  (fluxes) and registered in `Init()`. The fields staged to/from Noah-MP each
  step are listed in `NoahmpInputComp` (ERF → Noah-MP) and `NoahmpOutputComp`
  (Noah-MP → ERF). These are four *separate* lists (C++ `enum`s); keep them
  straight. See [`spec-noahmp-api.md`](spec-noahmp-api.md).
- **`noahmp_real` / `amrex::Real`.** ERF builds Noah-MP in double precision, and
  a run-time check (on the submodule side) aborts the run if the C++ and Fortran
  sides were compiled with different precision. The coupling precision follows
  `amrex::Real`, so never hardcode `double` in coupling code.

## Using these specs with a coding agent

These `dev/` documents are written to be loaded directly into a coding agent's
context (Claude Code or similar), not just read by people. The intended workflow:

1. **Give the agent this README first.** It is the map: which file owns which
   concern (the tables above) and which spec to open next.
2. **Then load the spec for the file being changed.** Each spec is scoped to a
   concern and cross-links the others, so the agent pulls in only what the task
   needs — API/lifecycle ([`spec-noahmp-api.md`](spec-noahmp-api.md)), per-step
   data movement ([`spec-noahmp-gpu.md`](spec-noahmp-gpu.md)), I/O
   ([`spec-noahmp-io.md`](spec-noahmp-io.md)), or source layout and the X-macro
   field registry ([`spec-noahmp-reorg.md`](spec-noahmp-reorg.md)).
3. **Treat the invariants as acceptance criteria.** The specs state what a change
   must not break; the agent should verify against them before finishing. For
   adding a coupled variable, the "Validation checklist" in
   [`spec-noahmp-api.md`](spec-noahmp-api.md) §"Adding a coupled variable" is the
   contract to check against.

This pattern — keeping design specs in `dev/` next to the code so humans and
agents work from the same contract — is part of ongoing AI-for-HPC research using
[CodeScribe](https://github.com/akashdhruv/CodeScribe).
