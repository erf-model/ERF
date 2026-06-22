# ERF Noah-MP driver — developer specs (`dev/`)

This directory holds the design specifications for the **ERF-side driver of
Noah-MP**: the `NOAHMP` C++ class in `ERF_NOAHMP.H` / `ERF_NOAHMP.cpp` that
plugs the Noah-MP land-surface model into ERF's land-surface-model interface and
drives the per-step state exchange between the two.

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

| File | What it covers |
|------|----------------|
| [`spec-noahmp-api.md`](spec-noahmp-api.md) | The `NOAHMP` class, its `NullSurf` contract, the `LsmData`/`LsmFlux` coupling fields, the run lifecycle, time subcycling, and the workflow for exposing a new coupled variable. **Start here.** |
| [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md) | The GPU-aware, component-indexed state exchange: pinned staging buffers, the `ParallelFor` (device) / `LoopOnCpu` (host) dataflow, synchronization points, the slab/ghost-cell handling, and the flux fill-value guard. |
| [`spec-noahmp-io.md`](spec-noahmp-io.md) | I/O owned by the driver: checkpoint/restart (the substep counter and the full prognostic state), the per-step land plotfile, and the static-input reads (`namelist.erf`, `NoahmpTable.TBL`, the NetCDF land file). |

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
- **CodeScribe.** `ERF_NOAHMP.cpp`/`.H` may optionally be updated with the
  LLM-based CodeScribe tool. The interface contract that prompt encodes is now
  captured in [`spec-noahmp-api.md`](spec-noahmp-api.md) §"Adding a coupled
  variable".
