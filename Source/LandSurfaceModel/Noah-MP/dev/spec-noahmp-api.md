# Spec: the `NOAHMP` driver class — physics coupling & API

> Status: living document · Owns: the `NOAHMP` C++ class, the ERF-side coupling
> contract, and the lifecycle. · Companion:
> [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md) (the per-step data movement),
> [`spec-noahmp-io.md`](spec-noahmp-io.md) (restart & output). · Files:
> `ERF_NOAHMP.H`, `ERF_NOAHMP.cpp`.

## 1. What this layer is

`NOAHMP` is ERF's concrete land-surface model for Noah-MP. It derives from
`NullSurf` (`Source/LandSurfaceModel/Null/ERF_NullSurf.H`), the abstract LSM
interface ERF programs against, and it owns one Fortran-backed `NoahmpIO_type`
per local box (the `noahmpio_vect` member). Everything ERF needs from a land
model — surface temperature and albedo for radiation, turbulent fluxes for the
surface layer, checkpoint/restart — is delivered through the `NullSurf` virtuals
this class overrides.

The driver does **not** contain land physics. It is the adapter between ERF's
AMReX data and the Noah-MP `NoahmpIO_type` API (documented in the submodule's
`drivers/erf/dev/` specs). Its responsibilities are: build the coupling
`MultiFab`s, drive Noah-MP on the right schedule, and shuttle state across the
ERF ↔ Noah-MP boundary each step (the shuttle itself is
[`spec-noahmp-gpu.md`](spec-noahmp-gpu.md)).

## 2. The four field enums (keep them straight)

Two *pairs* of enums describe two *different* boundaries. Confusing them is the
most common coupling bug.

### 2a. ERF-internal coupling fields (`LsmData_NOAHMP`, `LsmFlux_NOAHMP`)

These are the fields ERF's other components read from / write to the LSM, via the
`Lsm_Data_Ptr` / `Lsm_Flux_Ptr` accessors. Each is a 2-D (surface-slab)
`MultiFab` allocated in `Init()`.

- **`LsmData_NOAHMP`** — state shared with **radiation (RRTMGP)**: surface
  temperature `t_sfc`, emissivity `sfc_emis`, the four albedos
  (`sfc_alb_{dir,dif}_{vis,nir}`), `cos_zenith_angle`, and the downwelling
  short/long-wave fluxes that radiation produces and Noah-MP consumes
  (`sw_flux_dn*`, `lw_flux_dn`). Stored in `lsm_fab_data[]`.
- **`LsmFlux_NOAHMP`** — turbulent fluxes handed to the **surface layer**:
  `t_flux`, `q_flux`, and the momentum stresses `tau13`/`tau23`. Stored in
  `lsm_fab_flux[]` with one ghost cell — a one-wide halo of copied neighbor
  values — so the surface layer can average the stresses, which are staggered onto
  cell faces rather than cell centers.

Direction matters: some `LsmData` entries flow *into* Noah-MP as forcing
(`sw_flux_dn`, `lw_flux_dn`, `cos_zenith_angle`), the rest flow *out* of Noah-MP
as results (`t_sfc`, `sfc_emis`, albedos). Radiation runs first each step and
reads the *result* fields before Noah-MP has produced them, so those are
pre-seeded to sane values in `Init()` (`t_sfc = 300 K`, `sfc_emis = 0.9`,
albedos `= 0.06`); the remaining fields, including the forcing fields, are left
at zero.

Each enum has a `NumVars` sentinel. `Init()` builds parallel
`LsmDataMap`/`LsmDataName` (and the flux equivalents) so ERF can look fields up
by index or by name; the `*Map` arrays let ERF's generic index space map onto
this class's enum order.

### 2b. Per-step staging components (`NoahmpInputComp`, `NoahmpOutputComp`)

These index the *components* of the two multi-component pinned `FArrayBox`es used
to move data to/from Noah-MP each step. They are an implementation detail of the
exchange, **not** ERF-visible fields:

- **`NoahmpInputComp`** — forcing ERF → Noah-MP: `u_phy`, `v_phy`, `t_phy`,
  `qv_curr`, `p8w`, `swdown`, `glw`, `coszen`.
- **`NoahmpOutputComp`** — results Noah-MP → ERF: `hfx`, `lh`, `tau_ew`,
  `tau_ns`, `tsk`, `emiss`, and the four albedos.

Both end in `NumComps`, which sizes the staging buffers. The full mechanics of
how these components are written and read live in
[`spec-noahmp-gpu.md`](spec-noahmp-gpu.md).

## 3. The `NullSurf` contract this class implements

| Virtual | Role in `NOAHMP` |
|---------|------------------|
| `Define(sc)` | Hook for pulling constants from `SolverChoice` (currently a no-op). |
| `Init(lev, cons_in, geom, dt)` | Build coupling `MultiFab`s and `lsm` geometry, size + initialize one `NoahmpIO_type` per box, broadcast `DTBL`/counter. §4. |
| `Advance_With_State(...)` | Per-step forcing → `DriverMain()` → flux read-back, gated by the firing schedule. §5, and [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md). |
| `Lsm_Data_Ptr` / `Lsm_Flux_Ptr` | Return the coupling `MultiFab` for an ERF var index (via `LsmDataMap`/`LsmFluxMap`). |
| `Lsm_Data{Name,Index}` / `Lsm_Flux{Name,Index}` | Name ↔ index lookup over the enums. |
| `Lsm_Geom` / `Lsm_{Data,Flux}_Size` | Expose the surface geometry and field counts. |
| `Get_LSM_Step` / `Set_LSM_Step` | Read/write the substep counter for checkpoint/restart. §5, [`spec-noahmp-io.md`](spec-noahmp-io.md). |
| `Write_Lsm_Restart` / `Read_Lsm_Restart` | Serialize/restore the full Noah-MP prognostic state. [`spec-noahmp-io.md`](spec-noahmp-io.md). |
| `Plot_Landfile(nstep)` | Emit the per-step NetCDF land output. [`spec-noahmp-io.md`](spec-noahmp-io.md). |

## 4. Lifecycle: `Init()`

`Init()` runs once per AMR level and:

1. **Installs the Noah-MP fatal handler once** (across all levels/ranks). Noah-MP
   has no MPI of its own; routing `noahmp::fatal()` through `amrex::Abort`
   ensures a fatal on any rank becomes an `MPI_Abort` instead of deadlocking
   peers at the next collective. See `NoahmpFatal.H`.
2. **Builds the surface (`lsm`) geometry** — a single-cell-thick slab below the
   atmosphere's lowest level (`khi_lsm = domain.smallEnd(2) - 1`), with the same
   x/y decomposition as `cons_in` (every box ranged to `k = 0`).
3. **Registers the coupling fields** — fills `LsmDataMap`/`LsmDataName` and the
   flux equivalents, allocates `lsm_fab_data[]` / `lsm_fab_flux[]`, and seeds the
   result fields radiation reads before Noah-MP has produced them (§2a).
4. **Sizes the per-box state** — `noahmpio_vect.resize(local_size, lev)` (a rank
   owning no boxes leaves it empty), and allocates the pinned staging buffers
   `noahmp_input_tmp[]` / `noahmp_output_tmp[]` (one per box; see
   [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md)).
5. **Initializes each `NoahmpIO_type`** in the boundary-box loop:
   `ScalarInitDefault()`, `ReadNamelist()`, `ReadLandHeader()`, set the
   domain/memory/tile bounds, `VarInitDefault()`, `ReadTable()`,
   `ReadLandMain()`, `InitMain()`, then `WriteLand(0)`. (Method semantics belong
   to the submodule's API spec.)
6. **Broadcasts the firing parameters** so the schedule is well-defined on every
   rank, including land-free ranks. `m_dtbl` (Noah-MP's `DTBL`) and `m_itimestep`
   are reduced with `ReduceRealMax`/`ReduceIntMax`; land-free ranks seed
   sentinels (`lowest()`) that lose the max. Asserts `m_dtbl > 0` and
   `m_dt <= m_dtbl`.

> **Invariant — size once.** `noahmpio_vect` is sized exactly once. Each
> `NoahmpIO_type` object is self-referential — the Fortran side holds pointers
> back to the C++ scalars — so moving one in memory (which any later resize would
> do) leaves those pointers aimed at freed memory. Do not `push_back`/resize it
> later.

## 5. Firing schedule & time subcycling

The Noah-MP timestep (`NOAH_TIMESTEP` in `namelist.erf`, surfaced as `DTBL`) may
be larger than the ERF timestep, so Noah-MP subcycles: it fires once every
`DTBL/dt` ERF steps. The decision must be **identical on every MPI rank**, because
`Advance_With_State` ends in a `FillBoundary` ghost-cell exchange that every rank
must enter together (a *collective* operation). If one rank fired and another did
not, the run would deadlock at that exchange.

That is why the counter and `DTBL` are *class members broadcast to all ranks*
(`m_itimestep`, `m_dtbl`) rather than read from `noahmpio_vect[0]` (absent on
land-free ranks). Each call:

```
NOAH_time = (m_itimestep - 1) * m_dtbl;
if (elapsed_time < NOAH_time) return;   // not yet — every rank agrees
m_itimestep += 1;                        // advance in lockstep
```

The per-box `noahmpio->itimestep` is mirrored from `m_itimestep` just before
`DriverMain()`. `Get_LSM_Step`/`Set_LSM_Step` read/write `m_itimestep` (and fan
out to every box) so a restart resumes the schedule exactly — see
[`spec-noahmp-io.md`](spec-noahmp-io.md).

## 6. Adding a coupled variable

To thread a new field through the coupling, follow the memory flow used for
existing fields (`TSK`, `EMISS`, `HFX`, …). Preserve the loop structure, the
component-indexed buffer layout, the CPU↔GPU dataflow, and the existing variable
naming. The mechanics differ by direction:

**A. Forcing field (ERF → Noah-MP):**

1. Add an entry to `NoahmpInputComp` (`ERF_NOAHMP.H`) before `NumComps`. The
   buffers size off `NumComps`, so no buffer edit is needed.
2. In the device `ParallelFor` over `bx`, write the ERF source expression into
   that input component.
3. After `Gpu::streamSynchronize()`, in the host `LoopOnCpu` over `bx`, copy the
   component into the matching `noahmpio->FIELD(...)`, using the correct index
   rank (3-D atmospheric `(i,1,j)`; 2-D surface `(i,j)`; banded radiation
   `(i,band,j)`).

**B. Result field (Noah-MP → ERF):**

1. Add an entry to `NoahmpOutputComp` before `NumComps`.
2. In the post-`DriverMain()` host `LoopOnCpu`, read `noahmpio->FIELD(...)` into
   that output component.
3. In the device `ParallelFor` over `gbx`, assign the component into the
   destination ERF field, using the clamped `(ii,jj)` indices. A surface-layer
   flux divided by a state factor must follow the `t_flux`/`q_flux`/`tau13`/
   `tau23` pattern, including the `-9999` fill guard (see
   [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md) §"fill-value guard").

**C. New ERF destination field.** If the destination is not an existing
`lsm_fab_data` / `lsm_fab_flux` field, also add it to `LsmData_NOAHMP` /
`LsmFlux_NOAHMP` and register it in `Init()` (extend the `*Map`/`*Name` arrays
and the allocation/initialization loop).

**D. New Noah-MP-side member.** Exposing a brand-new `NoahmpIO_type` member (one
that does not exist on the Fortran side yet) is a *submodule* change — see the
`add-coupled-variable` spec under `Submodules/Noah-MP/drivers/erf/dev/`. The
steps above assume the `noahmpio->FIELD` member already exists.

**Validation checklist:**

- Array ranks/index ordering match (3-D → `(i,1,j)`, 2-D → `(i,j)`, banded →
  `(i,band,j)`; ERF surface `Array4` → `(i,j,0)`).
- Every new field has an enum entry **and** a matching read/write in each
  relevant copy block, so component indices stay one-to-one with `NumComps`.
- Sync order preserved: `Gpu::streamSynchronize()` before host access; the
  copy-back `ParallelFor` after `DriverMain()`.
- `noahmpio_vect` indexing (`idb`), the `bx.smallEnd(2) != klo` slab guard, and
  the `(ii,jj)` clamping are untouched.

> `ERF_NOAHMP.cpp`/`.H` can optionally be updated with **CodeScribe**
> (<https://github.com/akashdhruv/CodeScribe>), an LLM-based code-update tool.
> This section is the interface contract such a prompt must respect.
