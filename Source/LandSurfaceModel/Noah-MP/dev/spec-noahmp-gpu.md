# Spec: GPU-aware state exchange

> Status: living document · Owns: the per-step movement of data across the
> ERF (GPU) ↔ Noah-MP (host) boundary inside `Advance_With_State`. · Companion:
> [`spec-noahmp-api.md`](spec-noahmp-api.md) (the field enums & lifecycle),
> [`spec-noahmp-io.md`](spec-noahmp-io.md). · File: `ERF_NOAHMP.cpp`.

## 1. The problem

ERF state lives in device (GPU) memory and is updated by GPU kernels (AMReX
`ParallelFor` loops). Noah-MP is a host-only Fortran library: its `NoahmpIO_type`
arrays live in host (CPU) memory and `DriverMain()` runs on the CPU. Every step
the driver must move forcing fields *device → host*, run Noah-MP, then move
results *host → device* — without allocating a fresh buffer per variable and
without stalling the GPU more than necessary.

The solution is a pair of **pinned, multi-component staging buffers per box**.
"Pinned" host memory is page-locked so the GPU can write to it directly and the
CPU can read it; a single such buffer therefore bridges the two worlds.

## 2. The staging buffers

Allocated once per box in `Init()` and reused every step:

```cpp
noahmp_input_tmp[idb]  = std::make_unique<FArrayBox>(bx, NoahmpInputComp::NumComps , The_Pinned_Arena());
noahmp_output_tmp[idb] = std::make_unique<FArrayBox>(bx, NoahmpOutputComp::NumComps, The_Pinned_Arena());
```

Each is a single multi-component `FArrayBox` over the surface slab `bx`: one
coupled field per component, indexed by `NoahmpInputComp` / `NoahmpOutputComp`
(see [`spec-noahmp-api.md`](spec-noahmp-api.md) §2b). There is deliberately **no
per-variable `FArrayBox`** — adding a field means adding a component, not an
allocation, and the buffers size themselves off `NumComps`.

## 3. The six-step dataflow (per box, per firing)

`Advance_With_State` loops over boxes (`idb`) and, for each surface box, runs:

```
 (1) device ParallelFor over bx  ── write ERF forcing into input components ──┐
                                                                              │  pinned
 (2) Gpu::streamSynchronize()    ── device writes are visible to the host ────┤  input
                                                                              │  buffer
 (3) host LoopOnCpu over bx      ── copy input components → noahmpio->FIELD ──┘

 (4) noahmpio->itimestep = m_itimestep;  noahmpio->DriverMain();   ── run physics (host)

 (5) host LoopOnCpu over bx      ── copy noahmpio->FIELD → output components ──┐  pinned
                                                                              │  output
 (6) device ParallelFor over gbx ── copy output components → ERF fields ───────┘  buffer
```

Concretely:

- **(1) Device → input buffer.** A `ParallelFor` over `bx` computes each forcing
  field from ERF state and stores it in its input component, e.g.
  `noah_input_arr(i,j,0,NoahmpInputComp::t_phy) = getTgivenRandRTh(...)`. Wind is
  averaged from faces to cell centers; `qv`, `t_phy`, `p8w` are derived via the
  EOS helpers; `swdown`/`glw`/`coszen` come from the `lsm_fab_data` coupling
  fields.
- **(2) Synchronize.** `Gpu::streamSynchronize()` — **mandatory** before any host
  read of the pinned buffer; the device kernel is asynchronous.
- **(3) Input buffer → Noah-MP.** A `LoopOnCpu` over `bx` copies each component
  into the matching `noahmpio` member, applying the **k-axis transpose**:
  3-D atmospheric → `noahmpio->T_PHY(i,1,j)`; 2-D surface → `noahmpio->SWDOWN(i,j)`.
- **(4) Drive.** Mirror the authoritative counter (`noahmpio->itimestep =
  m_itimestep`) then `DriverMain()`.
- **(5) Noah-MP → output buffer.** A `LoopOnCpu` over `bx` reads each result into
  its output component. Banded radiation reads use the band index:
  `ALBSFCDIRXY(i,1,j)` (VIS), `ALBSFCDIRXY(i,2,j)` (NIR).
- **(6) Output buffer → device.** A `ParallelFor` over the grown box `gbx` writes
  the results into the destination ERF fields (`lsm_fab_data` for RRTMGP vars,
  `lsm_fab_flux` for surface-layer fluxes), with index clamping and the
  fill-value guard (§5).

## 4. Index conventions

| Layout | Addressing | Example |
|--------|------------|---------|
| ERF `Array4` (state, staging buffers, coupling fields) | `(i, j, k)` / `(i, j, 0)` on the slab | `CONS(i,j,k,Rho_comp)`, `noah_input_arr(i,j,0,comp)`, `TSK(i,j,0)` |
| Noah-MP 3-D atmospheric | `(i, k, j)`, k pinned to surface layer `1` | `noahmpio->T_PHY(i,1,j)` |
| Noah-MP 2-D surface | `(i, j)` | `noahmpio->TSK(i,j)` |
| Noah-MP banded radiation | `(i, band, j)`, VIS=1 NIR=2 | `noahmpio->ALBSFCDIRXY(i,1,j)` |

The **k/j transpose** between ERF `(i,j,k)` and Noah-MP `(i,k,j)` is the single
most error-prone part of the exchange. Getting it wrong does not crash; it
silently corrupts the column.

## 5. Slab guard, ghost cells, and the flux fill-value guard

- **Surface-slab guard.** Only the box touching the domain floor participates:
  `if (bx.smallEnd(2) != klo) continue;` then `bx.makeSlab(2, klo)`. Boxes higher
  in the column are skipped. This guard appears in both `Init()` and
  `Advance_With_State` and must stay in sync.
- **Grown box for copy-back.** Stage (6) runs over `gbx` (grown by one in x/y) so
  ghost cells get values; indices are clamped to the valid box
  (`ii = clamp(i, i_lo, i_hi)`, likewise `jj`) so the ghost ring copies the
  nearest interior result. After the box loop, `lsm_fab_flux[*]->FillBoundary`
  reconciles the ghost cells across boxes (`tau13`/`tau23` are staggered onto cell
  faces, so the surface layer needs the one-cell halo to average them back to cell
  centers).
- **Fill-value guard (fluxes only).** Noah-MP returns `-9999` for cells it does
  not process (sea-ice / open-water points that still carry `LANDMASK=1`).
  Applying that as a flux (`-9999/(rho*Cp)`) would crash the lowest cell to
  ~200 K. So stage (6) tests `hfx_lsm > -9990` before dividing; failing cells get
  the `lsm_flux_undefined` sentinel, and the surface layer falls back to its MOST
  flux there (see `ERF_SurfaceLayer.cpp`). The four fluxes
  (`t_flux`/`q_flux`/`tau13`/`tau23`) are written together under this guard.

## 6. Invariants (do not regress)

1. **Sync before host read.** `Gpu::streamSynchronize()` must sit between the
   device write (stage 1) and the host read (stage 3). Removing it is a data
   race that only shows up under load.
2. **One component per field, one-to-one with `NumComps`.** Every field has an
   enum entry *and* a read/write in each relevant copy block. A gap desyncs all
   higher component indices.
3. **Buffers sized off `NumComps`, never hand-sized.** Adding a field is an enum
   line plus copy-block lines; never edit the `make_unique<FArrayBox>` extent.
4. **Transpose respected.** 3-D Noah-MP access is `(i,1,j)`, not `(i,j,1)`.
5. **Pinned arena.** Staging buffers come from `The_Pinned_Arena()`; ordinary
   device or managed memory breaks the host read in stage (3)/(5).
6. **Counter mirrored, not owned, per box.** `noahmpio->itimestep` is set from
   the broadcast `m_itimestep` each firing (see
   [`spec-noahmp-api.md`](spec-noahmp-api.md) §5); the per-box value is not the
   source of truth.
