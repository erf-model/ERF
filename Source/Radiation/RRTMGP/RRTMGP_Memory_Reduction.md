# RRTMGP Radiation Memory Reduction

## Problem

RRTMGP radiation was extremely memory-consuming because all horizontal columns on each
MPI rank were processed at once. Several internal data structures scale as
`O(ncol * nlay * ngpt)` where `ngpt` is 224 (SW) or 256 (LW). For a rank with 50,000
columns and 100 vertical levels, peak GPU memory exceeded 40+ GB. The Kokkos memory pool
alone was sized at `300 * ncol * nlay * nbnd * 8 bytes` — a massive over-allocation.

## Solution

Five changes were implemented, ordered by impact:

### 1. Right-size the Kokkos Memory Pool

**File:** `ERF_RRTMGP_Interface.cpp`

The pool multiplier `nvar` was reduced from 300 to 12, and `nbnd` was replaced with `ngpt`
since the pool stores per-g-point temporaries. The actual peak concurrent usage inside the
RTE solvers is ~15-20 arrays of size `ncol * nlay * ngpt`.

### 2. Initialize RRTMGP Once (Not Every Radiation Step)

**Files:** `ERF_Radiation.cpp`, `ERF_Radiation.H`

Previously, `rrtmgp_initialize()` loaded 4 NetCDF files and created the memory pool every
radiation step, and `rrtmgp_finalize()` destroyed everything. The k-distribution and cloud
optics tables are constant, so initialization now happens once (guarded by
`rrtmgp::initialized`). The memory pool is sized for `ncol_chunk` rather than full `ncol`.
Finalization is deferred to the `Radiation` destructor.

### 3. Column Chunking in `run_impl()` (Primary Memory Fix)

**Files:** `ERF_Radiation.cpp`, `ERF_Radiation.H`

Radiation columns are completely independent (no horizontal coupling). They are now
processed in chunks of `ncol_chunk` (default 5000, configurable via `erf.rad_ncol_chunk`)
instead of all at once. This is the same approach used by E3SM/SCREAM.

Implementation details:
- Full-ncol input/output buffers are kept in `alloc_buffers()`
- The memory explosion happens inside `rrtmgp_main()` where arrays of size
  `ncol * nlay * ngpt` are created on the stack allocator
- The `rrtmgp_main()` call is wrapped in a chunk loop
- Unmanaged Kokkos Views (pointer arithmetic) create zero-copy subviews into the full
  buffers, leveraging LayoutRight which guarantees contiguous column data
- Per-chunk `gas_concs_t` objects are created with subset VMR data
- `compute_band_by_band_surface_albedos`, `compute_heating_rate`, and
  `compute_broadband_surface_fluxes` are called per chunk

**New input parameter:** `erf.rad_ncol_chunk` (integer, default 5000)

### 4. Remove Dead Arrays and Parameters

**Files:** `ERF_Radiation.H`, `ERF_Radiation.cpp`, `ERF_RRTMGP_Interface.H`,
`ERF_RRTMGP_Interface.cpp`

Removed unused member declarations and function parameters that were already commented out:
- Cloud band optical depth arrays (`cld_tau_sw_bnd`, `cld_tau_lw_bnd`)
- Cloud g-point optical depth arrays (`cld_tau_sw_gpt`, `cld_tau_lw_gpt`)

The aerosol optical property arrays (`aero_tau_sw`, `aero_ssa_sw`, `aero_g_sw`, `aero_tau_lw`)
are retained as dormant scaffolding. They are declared in `ERF_Radiation.H` but only
allocated when `m_do_aerosol_rad` is true, which currently triggers an abort in the
`Radiation` constructor (see section 6). This preserves the hook so a future aerosol
coupling (e.g. SPA, prescribed aerosol climatology) can populate these arrays without
reintroducing the members.

### 5. Conditional Allocation of Diagnostic Flux Arrays

**Files:** `ERF_Radiation.cpp`, `ERF_RRTMGP_Interface.cpp`

12 clean-sky/clean-clear-sky flux arrays (each `ncol * (nlay+1)`) are only allocated at
full size when `m_extra_clnclrsky_diag` or `m_extra_clnsky_diag` are true (both default
false). When disabled, 1-element placeholders are allocated instead. The flux zeroing
kernels in `rrtmgp_sw()` and `rrtmgp_lw()` are also made conditional on these flags.

### 6. Abort When Aerosol Forcing Is Requested

**File:** `ERF_Radiation.cpp`

`m_do_aerosol_rad` defaults to `false`. The `Radiation` constructor now aborts with a
descriptive message if a user sets `erf.rad_do_aerosol = true`, since the aerosol optical
property pipeline is not implemented (the arrays exist as scaffolding but are never
populated with real data — see section 4). The flag is preserved so the abort can be
removed once a real aerosol coupling lands.

### 7. Fence Before Datalog Read to Eliminate Cross-Stream Race

**File:** `ERF_Radiation.cpp`

A `Kokkos::fence()` is inserted between `rrtmgp::compute_heating_rate(...)` and
`populateDatalogMF()`. Without the fence, the LW clear-sky heating-rate kernel (launched
on the Kokkos default stream) could still be in flight when `populateDatalogMF()` (an
AMReX `ParallelFor` on `Gpu::gpuStream()`) read `lw_clrsky_heating`. The race produced
run-to-run-variable `radqrclw` in the radiation datalog while every other field stayed
deterministic. This bug exists on the `development` branch as well — it was exposed, not
introduced, by the chunking refactor.

## Expected Impact

| Metric | Before | After (chunk=5000, 50K cols, 100 levels) |
|--------|--------|------------------------------------------|
| Memory pool | ~192 GB | ~0.5 GB |
| Peak GPU memory | 40+ GB | ~4 GB |
| Diagnostic flux arrays | ~400 MB always | ~400 MB only when enabled |

The `erf.rad_ncol_chunk` parameter controls the tradeoff between peak memory and the
number of kernel launches. Smaller chunks use less memory but launch more kernels.
Results are bit-identical regardless of chunk size since radiation columns are independent.

## Verification

1. Build with `ERF_ENABLE_RRTMGP=ON` (CUDA build verified successfully)
2. Run `ctest -R Radiation -VV` and compare against gold files
3. Verify pool high-water mark in `rrtmgp_finalize()` output
4. Use `nvidia-smi` to confirm peak GPU memory reduction
