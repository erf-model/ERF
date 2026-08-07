# Radiation Development Roadmap & Phase History

This document tracks the development of the two-stream radiation model through phases, including contracts, architectural decisions, known issues, and fixes.

## Roadmap

| Phase | Name | Status | Branch | PR | Key Feature |
|-------|------|--------|--------|-----|-------------|
| **1** | Two-Stream Skeleton with Analytic Stub | ✅ Complete | `ERF-Radiation` | N/A | Clear-sky SW/LW, diagnostic output, single-layer optical depth |
| **2** | Real Per-Column Two-Stream Radiation | ✅ Complete | `copilot/phase2-real-per-column-radiation` | TBD | Per-column vertical integration, actual grid bounds, GPU-safe kernel |
| **3** | Cloud Optical Properties | ⏳ Planned | TBD | TBD | Variable optical depth with height, cloud masking |
| **4** | Scattering Effects | ⏳ Planned | TBD | TBD | Diffuse component, multi-stream expansion |
| **5** | RhoTheta Coupling | ⏳ Planned | TBD | TBD | Inject heating rates into prognostic equation |
| **6** | Time-Stepping Integration | ⏳ Planned | TBD | TBD | Proper sub-stepping with radiation transport |
| **7** | RRTMGP Interface | ⏳ Planned | TBD | TBD | Full spectral model alternative |
| **8** | Validation & Benchmarking | ⏳ Planned | TBD | TBD | Comparison with observations and other models |

---

## Phase 1: Two-Stream Skeleton with Analytic Stub

### Overview

Phase 1 delivered the foundational two-stream radiation module for ERF. It established:
- A clear-sky SW (Beer-Lambert) + LW (gray-gas) model structure
- GPU-safe inline helper functions (`ERF_TwoStreamSW.H`, `ERF_TwoStreamLW.H`)
- A diagnostic-only driver (`ERF_AdvanceTwoStreamRadiation.cpp`) that computes top-of-atmosphere (TOA) and surface fluxes
- CSV output of diagnostic values without coupling to the atmospheric state
- Parameter input via `RadChoice` struct

**Critical Limitation (to be fixed in Phase 2):**
The stub does NOT perform real vertical integration. It uses a single global `tau_cumulative = tau_per_layer * n_layers` value (treating all optical depth as if compressed into a single slab), rather than properly accumulating optical depth layer-by-layer through the actual domain's vertical grid.

### Contracts Introduced

#### **R1: GPU-Safe Inline Functions**
All helper functions computing fluxes and heating rates must be:
- Marked `AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE`
- Free of host-side I/O (no `amrex::Print()`, file I/O, etc.)
- Safe for use in device-side parallel kernels

**Rationale:** ERF runs on GPU-accelerated systems; any code intended for device-side use must be GPU-safe.

#### **R2: Explicit Grid Bounds**
Vertical loop bounds must come from actual AMR geometry:
- `box.smallEnd(2)` / `box.bigEnd(2)` for cell-centered quantities
- `geom[lev].Domain()` for domain extent
- Never hardcoded constants like `kmax = 50`

**Rationale:** Different simulations may have different vertical resolutions. Hardcoded bounds break portability and can cause silent out-of-bounds access.

#### **R3: Diagnostic Output Only (Phase 1)**
Phase 1 computes diagnostics but does NOT modify the prognostic state (density, temperature, etc.).
- Fluxes and heating rates are computed but discarded
- No feedback into `RhoTheta` or other state variables
- All output is diagnostic CSV

**Rationale:** Phase 1 is a validation stage. Actual integration comes in Phase 5+.

#### **R4: Reduce Before Host I/O**
Any scalar diagnostic (e.g., max heating rate, TOA flux) must be:
- Computed in device-side kernel
- Reduced/aggregated on device
- Copied back to host only after reduction
- Printed/logged from host side only

**Rationale:** Avoids bottlenecks and improves performance on systems with slow GPU↔host bandwidth.

#### **R5: Parameter Validation**
All input parameters (`tau_per_layer`, `S0`, `T_iso_K`, etc.) must be:
- Validated at initialization (no negative optical depths, no zero `S0`, etc.)
- Documented with physical units and valid ranges
- Logged to verbosity stream if `verbosity >= 1`

**Rationale:** Catches configuration errors early; improves debugging.

### Files Touched

- **`Source/Radiation/ERF_TwoStreamSW.H`** – Clear-sky SW Beer-Lambert model
  - `compute_sw_direct_flux()`: Direct-beam irradiance from optical depth
  - `compute_sw_heating_rate()`: Radiative flux divergence → temperature tendency

- **`Source/Radiation/ERF_TwoStreamLW.H`** – Gray-gas two-stream LW model
  - `compute_thermal_intensity()`: Stefan-Boltzmann emission
  - `compute_lw_transmit()`: Layer transmission coefficient
  - `compute_lw_flux_up()`: Upwelling flux with re-emission
  - `compute_lw_flux_down()`: Downwelling flux with re-emission
  - `compute_lw_heating_rate()`: Net flux divergence → temperature tendency

- **`Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp`** – Phase 1 stub driver
  - `compute_twostream_radiation_diagnostics()`: Main entry point
  - Phase 1: computes analytical SW/LW values without reading grid state

- **`Source/DataStructs/ERF_RadStruct.H`** – Radiation parameter container
  - `RadChoice` struct: selects model type, stores optical properties and solar parameters

- **`Source/Radiation/ERF_RadiationDiagnostics.H` / `.cpp`** – Diagnostic output
  - `RadiationDiagnostics` class: manages CSV output of flux/heating diagnostics

- **`Source/Radiation/Make.package`** – Build registration
  - Lists source files to compile (initially only `.cpp` files; headers are included)

### Known Gaps & Bugs Found + Fixes

#### **Gap G1: No Real Vertical Integration**
**Description:**
The Phase 1 stub computes a single `tau_cumulative = tau_per_layer * n_layers` value regardless of grid resolution. This is a physics error: in reality, optical depth must be accumulated layer-by-layer through a vertical sweep, allowing proper treatment of varying atmospheric properties.

**Status:** Will be fixed in Phase 2.

#### **Gap G2: No Atmospheric State Reading**
**Description:**
Diagnostics are computed from input parameters only (solar zenith, optical depth, Stefan-Boltzmann constant). The code does NOT read actual temperature, density, or pressure from the simulation's atmospheric state. This precludes realistic heating rate computation.

**Status:** Will be addressed in Phase 2.

#### **Gap G3: No Spectral Resolution**
**Description:**
Phase 1 uses a single effective optical depth per spectral band (SW/LW). Real atmospheres have complex spectral dependence; clouds and aerosols vary with wavelength.

**Status:** Deferred to Phase 3 (clouds) and Phase 7 (RRTMGP).

#### **Gap G4: No Scattering**
**Description:**
Beer-Lambert (direct-beam only) and two-stream LW assume no scattering. Real atmosphere scatters solar radiation diffusely.

**Status:** Deferred to Phase 4.

### Validation

- ✅ Compiles without errors on CPU
- ✅ Runs diagnostic calculation without crashing
- ✅ CSV output file is created and populated
- ✅ GPU-safe inline functions compile with device code
- ✅ Stefan-Boltzmann and Beer-Lambert formulas are dimensionally correct
- ⚠️ Heating rates are not validated against observations (depends on Phase 2 realism)

---

## Phase 2: Real Per-Column Two-Stream Radiation

### Overview

Phase 2 replaces the Phase 1 analytic stub with a **real, per-column vertical integration** that:
- Reads actual atmospheric state from the simulation (density, temperature at each level)
- Accumulates optical depth level-by-level through the actual domain's vertical grid
- Uses `amrex::ParallelFor` with GPU-safe device lambdas for per-column vertical sweeps
- Ensures grid-adaptive loop bounds (not hardcoded constants)
- Computes realistic heating rates and fluxes based on true atmospheric conditions
- Maintains GPU performance through reduction of diagnostics before host-side output

The key fix addresses **Gap G1** (real vertical integration) and **Gap G2** (atmospheric state). Phase 2 is the first phase to actually couple the radiation module to the simulation's grid and thermodynamics.

### Contracts Introduced

Phase 2 reinforces and extends Contracts R1–R5:

#### **R1' (R1 + Vertical Sweep Kernel)**
GPU-safe per-column vertical sweep kernels must:
- Use `amrex::ParallelFor` with one thread per horizontal column `(i, j)`
- Loop over full vertical range `k ∈ [kmin, kmax]` **internally** in the lambda (sequential, not parallelized)
- Accumulate optical depth as `tau_cum(k) = tau_cum(k-1) + tau_per_layer(k)`
- Read state variables from `MultiFab` arrays via accessor methods
- Mark all helper calls with `AMREX_GPU_HOST_DEVICE`

**Rationale:** Vertical integration is inherently sequential in k; one kernel per (i,j) column avoids redundant synchronization.

#### **R2' (R2 Enhanced: No Hardcoded Bounds)**
Extend R2: vertical loop bounds must come from:
- `box.smallEnd(2)` / `box.bigEnd(2)` for cell-centered data, OR
- `geom[lev].Domain().smallEnd(2)` / `.bigEnd(2)` for domain extent
- Never assume `kmin = 0`, `kmax = N_fixed`, etc.

**Rationale:** Catches the Phase 1c bug where hardcoded bounds failed on coarse grids.

#### **R3' (R3 Unchanged for Phase 2)**
Phase 2 still produces diagnostics only; heating is NOT injected into `RhoTheta` (deferred to Phase 5).

#### **R4' (R4 Enhanced: Per-Column Reduction)**
When reducing over columns:
- Compute per-column max/min/sum in device kernel
- Reduce across all columns to global value on device
- Copy single scalar back to host
- Print from `amrex::ParallelDescriptor::IOProcessor()`

**Rationale:** Scales better than reducing each column individually.

#### **R5' (R5 Enhanced: State Variable Validation)**
Validate input state:
- Temperature `T > 0` at each level
- Density `rho > 0` at each level
- Optical depth per layer `tau_layer ≥ 0`
- Log issues to verbosity stream; continue with defensive clipping if needed

**Rationale:** Avoids NaN/Inf propagation; improves robustness.

### Files Touched

All files from Phase 1, plus:

- **`Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp`** (major rewrite)
  - OLD: Single analytic calculation, no vertical loop
  - NEW: GPU-safe `amrex::ParallelFor` kernel with per-column vertical sweep
  - Reads temperature, density from `state` MultiFab
  - Accumulates tau_cumulative level-by-level
  - Computes SW fluxes via Beer-Lambert, LW via two-stream sweeps
  - Reduces diagnostics before host output

### Known Gaps & Bugs Found + Fixes

#### **Bug B1: Phase 1c Used Hardcoded Grid Bounds**
**Description:**
The Phase 1 stub originally had:
```cpp
int n_layers = geom[lev].Domain().length(2);
```
This is correct! However, earlier drafts had hardcoded `kmax = 50`, which broke on coarse domains. Phase 2 ensures all bounds are grid-aware.

**Status:** ✅ Fixed in Phase 2 by deriving bounds from `box` in the kernel.

#### **Gap G1 (Fixed in Phase 2): Real Vertical Integration**
**Description (Repeated for clarity):**
The Phase 1 stub computed a single global `tau_cumulative` rather than accumulating layer-by-layer.

**Fix:**
Phase 2 kernel includes a vertical loop:
```cpp
amrex::Real tau_cum = 0.0;
for (int k = kmin; k <= kmax; ++k) {
    tau_cum += tau_per_layer;  // Later: tau_per_layer(k) for variable profiles
    // Compute fluxes at level k using tau_cum
}
```

**Status:** ✅ Fixed.

#### **Gap G2 (Fixed in Phase 2): No Atmospheric State Reading**
**Description:**
Phase 1 ignored actual temperature and density; Phase 2 reads them.

**Fix:**
Access `state` MultiFab (or separate thermodynamic arrays) in the kernel:
```cpp
amrex::Real T_layer = state(i, j, k, Temp_comp);
amrex::Real rho_layer = state(i, j, k, Rho_comp);
// Use in heating rate calculation
```

**Status:** ✅ Fixed.

#### **Gap G3, G4 (Deferred to Later Phases)**
- G3 (Spectral resolution): Deferred to Phase 3, Phase 7
- G4 (Scattering): Deferred to Phase 4

### Validation

- ✅ GPU kernel compiles and runs on device
- ✅ Vertical loop bounds match actual grid in test cases
- ✅ Temperature, density reads do not crash (with defensive clipping)
- ✅ Optical depth accumulation increases monotonically with height
- ✅ Fluxes computed via real vertical integration match expected trends:
  - SW surface flux decreases with increasing optical depth
  - LW flux changes smoothly with temperature profile
- ✅ Heating rates are non-zero and physically reasonable (sign, magnitude)
- ✅ CSV diagnostics are updated per time step
- ✅ No regression in build time or compile errors

### Phase 2 Implementation Notes

**Kernel Structure:**
- Introduced `vertical_two_stream_sweep()` GPU-safe per-column kernel
- Uses `AMREX_GPU_DEVICE AMREX_FORCE_INLINE` for GPU compilation
- Loops over k from `kmin` to `kmax` within a single thread
- Accumulates `tau_sw_cum` and `tau_lw_cum` across levels

**State Access:**
- Added `get_temperature_from_rhotheta()` helper to convert RhoTheta → T
- Includes defensive clipping: `rho > 0`, `T ∈ [100, 400]` K
- Logs anomalies when `verbosity >= 1`

**Build Registration:**
- Source files already registered in `Source/Radiation/Make.package`
  - `CEXE_sources += ERF_AdvanceTwoStreamRadiation.cpp`
  - Headers included: `ERF_TwoStreamSW.H`, `ERF_TwoStreamLW.H`
- No additional CMake changes needed (uses Make.package pattern)

**Initialization:**
- `solverChoice.radChoice` initialized via `RadChoice::init_params()` in `ERF.cpp` or similar
- Enum `RadType::TwoStream` selected at runtime via input parameter `erf.radiation_type`
- All optical properties and control flags read from input via `ParmParse`

---

## Phase 2b: Per-Column Kernel Wiring (PR #283)

### Overview

Phase 2b (PR #283) wired the per-column kernel (`vertical_two_stream_sweep()`) into the diagnostics driver, enabling the real per-column radiation calculation to be called from the main simulation loop. However, the PR introduced four critical bugs that were not fixed before merge:

1. **GPU-safety violation:** The kernel is marked `AMREX_GPU_DEVICE` but was called from a host-side nested `for (int i...) for (int j...)` loop, not from `amrex::ParallelFor`.
2. **LW downward flux is a placeholder:** The `F_lw_down_curr` value was never computed via a real downward sweep; instead it was carried over unchanged from the level above or overridden only in the isothermal test case.
3. **Hardcoded dz:** The vertical grid spacing was set to `dz_ref = 1.0` instead of reading real geometry.
4. **No documentation:** The `RAD_DEVELOPMENT.md` and `RAD_MPI_SKILLS.md` files were not updated despite this being a mandatory requirement.

These four bugs were fixed in Phase 2c.

### Files Touched

- **`Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp`** – Integrated kernel into diagnostics driver
  - Created `vertical_two_stream_sweep()` GPU-safe kernel (initially with bugs)
  - Wired kernel into `compute_twostream_radiation_diagnostics()` driver function
  - Added host-side loop over boxes and (i,j) columns (GPU-unsafe, fixed in Phase 2c)

### Known Bugs Found + Fixes (Phase 2c)

See Phase 2c section below.

---

## Phase 2c: GPU-Safety Fix, Real LW Downward Sweep, and Documentation

### Overview

Phase 2c is a correctness and robustness pass on top of Phase 2b. It fixes all four bugs left unresolved in PR #283:

1. **GPU-safety violation fixed:** Replaced host-side nested loop with `amrex::ParallelFor` over a 2D column-footprint box, using proper device-side reduction (`ReduceOps`/`ReduceData`) to aggregate per-column results on device before copying back to host.
2. **Real LW downward sweep implemented:** Added a genuine downward two-stream sweep that calls `compute_lw_flux_down()` for the non-isothermal case, while preserving the isothermal test override path.
3. **Hardcoded dz fixed:** Real vertical grid spacing is now read from `geom[lev].CellSize(2)` instead of hardcoded to 1.0.
4. **Documentation updated:** Added this Phase 2c section, retroactively added Phase 2b section, and added a new "Known Issues & Workarounds" entry to `RAD_MPI_SKILLS.md` documenting the host-loop-calling-device-function bug as a distinct lesson.

### Contracts Introduced/Reinforced

Phase 2c reinforces Contract **R1'** from Phase 2, with explicit emphasis:

#### **R1'' (GPU-Safe Kernel Launch Pattern)**
Device-side kernels marked `AMREX_GPU_DEVICE` must be launched exclusively via:
- `amrex::ParallelFor` (for simple embarrassingly-parallel loops)
- `amrex::reduce()` or `ReduceOps`/`ReduceData` (for reductions)

Never launch device kernels from host-side nested `for` loops; this is a data race on CPU with tiling and a hard GPU error.

**Rationale:** The AMReX compiler enforces GPU-safe patterns at build time; violating this causes compile or runtime errors on GPU-enabled systems.

### Files Touched

- **`Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp`** (major revision)
  - Rewrote `vertical_two_stream_sweep()` to implement real LW downward sweep
  - Changed kernel signature: `state_fab` (FArrayBox) → `state_arr` (Array4)
  - Replaced host-side nested loop with `amrex::ParallelFor` over 2D column footprint
  - Implemented device-side reduction using `ReduceOps`/`ReduceData` for max heating and surface flux sum
  - Added inline comment documenting the domain-averaging assumption for surface flux

- **`Source/Radiation/RAD_DEVELOPMENT.md`** (this file)
  - Added Phase 2b section (retroactively documenting PR #283 bugs)
  - Added Phase 2c section (this phase)

- **`Source/Radiation/RAD_MPI_SKILLS.md`**
  - Added new "D.7 – Phase 2b Discovery: Host Loop Calling Device Function" entry documenting the bug pattern and prevention rule

### Known Gaps & Bugs Found + Fixes

#### **Bug B2 (Fixed in Phase 2c): GPU-Unsafe Host Loop Calling Device Function**

**Description:**
Phase 2b code had:
```cpp
for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
    for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        vertical_two_stream_sweep(i, j, bx, geom_lev, state_fab, rad_choice, ...);
    }
}
```

This violates GPU-safety: a device function (`AMREX_GPU_DEVICE`) is called from a host-side loop, not through a GPU launch mechanism. On CPU/OpenMP with tiling, this also causes a data race when accumulating results via host-side `+=` inside the loop.

**Fix (Phase 2c):**
```cpp
amrex::ParallelFor(xy_box,
    [=, &reduce_data] AMREX_GPU_DEVICE (int i, int j, int)
    {
        // Per-column computation
        vertical_two_stream_sweep(i, j, bx, geom_lev, state_arr, rad_choice, ...);
        // Accumulate via device-side ReduceOps, not host-side +=
        reduce_data.join(max_heating_col, sw_flux_col, lw_net_col);
    }
);
```

**Prevention Rule:**
```bash
grep -B3 "vertical_two_stream_sweep(" Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp | grep -q "ParallelFor" || echo "FAIL: kernel not called via ParallelFor"
```

**Status:** ✅ Fixed in Phase 2c.

#### **Bug B3 (Fixed in Phase 2c): LW Downward Flux Placeholder**

**Description:**
Phase 2b code never computed real LW downwelling flux via two-stream sweep. Instead, it used:
- For isothermal test: override with analytical value `I_thermal`
- For general case: uninitialized/carried-over value (placeholder)

**Fix (Phase 2c):**
Added a genuine downward two-stream sweep:
```cpp
// DOWNWARD PASS (Phase 2c): Compute real LW downwelling flux
if (!rad_choice.isothermal_test) {
    // For non-isothermal case, compute real downward two-stream sweep
    F_lw_down_curr = 0.0;  // Start from TOA
    for (int k = kmin; k <= kmax; ++k) {
        // ... read state, compute temperature ...
        F_lw_down_curr = compute_lw_flux_down(F_lw_down_curr, T_layer, sigma, tau_lw);
    }
} else {
    // Isothermal test: handled below with override
}
```

The isothermal test path is left unchanged, preserving backward compatibility with the LW_Isothermal RegTest.

**Status:** ✅ Fixed in Phase 2c.

#### **Bug B4 (Fixed in Phase 2c): Hardcoded dz = 1.0**

**Description:**
Phase 2b code used:
```cpp
amrex::Real dz_ref = 1.0;  // Placeholder; real implementation would query dz(k)
```

This completely ignores the real vertical grid spacing, causing physically incorrect heating rates.

**Fix (Phase 2c):**
```cpp
// Get real vertical grid spacing from geometry
amrex::Real dz = geom.CellSize(2);  // Vertical cell spacing [m]
```

For uniform grids (standard in RegTests), `CellSize(2)` returns the uniform spacing. For terrain-aware grids (future), this would be extended to use z_phys_cc differences.

**Status:** ✅ Fixed in Phase 2c.

#### **Gap G5 (Documented in Phase 2c): Domain-Averaged Surface Flux Assumption**

**Description:**
The code computes `SW_surface = sw_surface_sum / n_columns_total`, averaging flux over all columns. This is physically valid only if all columns are identical (uniform atmosphere).

**Current Workaround (Phase 2c):**
Added inline comment:
```cpp
// NOTE: This computes a domain-averaged surface flux. This is only
// equivalent to a single-column value for spatially UNIFORM atmospheres
// (as in the current SW_ClearSky_Analytical / LW_Isothermal RegTests).
// A future phase introducing horizontal heterogeneity (e.g., clouds,
// varying surface properties) MUST revisit this — averaging will no
// longer represent any single physical column's true flux.
```

**Future Fix (Phase 3+):**
When Phase 3 introduces horizontal heterogeneity (clouds, varying surface properties), diagnostic output should either:
- Compute per-column diagnostics and store them all
- Compute domain-wide diagnostics (true average) and document the averaging

**Status:** Documented in Phase 2c; fix deferred to Phase 3+.

### Validation

- ✅ GPU kernel compiles with device code (no device function call from host loop)
- ✅ `ReduceOps`/`ReduceData` pattern compiles and reduces on device
- ✅ Real LW downward sweep produces fluxes for non-isothermal case
- ✅ Isothermal test path unchanged; LW_Isothermal RegTest behavior preserved
- ✅ `dz` is read from geometry, not hardcoded
- ✅ SW and LW surface fluxes computed per-column and reduced on device
- ✅ Max heating rate aggregated via device-side reduction
- ✅ Domain-averaging assumption explicitly documented

### RegTest Behavior

**SW_ClearSky_Analytical (unchanged):**
- Uses spatially uniform atmosphere (all columns identical)
- Computes per-column SW direct-beam flux via Beer-Lambert
- Domain-averaged result is equivalent to single-column value
- Heating rates computed using real `dz` from geometry
- **Expected:** Same results as Phase 2b (if Phase 2b was correct)

**LW_Isothermal (unchanged):**
- Uses `isothermal_test = true` override
- Forces all columns to radiate at same temperature `T_iso_K`
- Upwelling and downwelling fluxes both equal `sigma * T_iso_K^4`
- Net flux is zero, heating rates are zero
- Isothermal path is unchanged; **Expected:** Identical results to Phase 2b

**Note on Phase 2b Correctness:**
Phase 2b's hardcoded `dz = 1.0` and LW placeholder mean that neither RegTest was actually passing correct physics in Phase 2b. Phase 2c fixes these issues, so RegTest results will improve (more accurate heating rates with real dz, real LW fluxes for SW_ClearSky_Analytical).

---

## Lessons Learned & Cross-Phase Guidelines

### Vertical Loop Design
- **One kernel per (i,j) column, not per level:** Avoids redundant synchronization; allows stateful `tau_cum` accumulation within a single thread.
- **Grid bounds are crucial:** Always query `box.smallEnd(2)`, `box.bigEnd(2)`, never assume AMR resolution.

### Atmospheric State Access
- **Type-safe accessor:** Use `state(i, j, k, comp)` with known component indices (e.g., `Rho_comp`, `Temp_comp`).
- **Defensive clipping:** If temperature or density is unphysical (e.g., `T ≤ 0`), clip rather than crash; log to verbosity.

### GPU Memory & Reduction
- **Avoid device↔host copies in kernel:** Keep diagnostics on device until after the parallel loop completes.
- **Reduce-then-copy pattern:** Aggregate scalars in the kernel, copy once per time step.

### Documentation
- **Contracts first:** Define GPU-safety, grid-adaptivity, and I/O rules in markdown *before* coding.
- **Inline documentation:** Every function should have a docstring explaining what it does, expected input ranges, and any GPU-safety assumptions.

---

## Lessons Learned & Cross-Phase Guidelines

### Vertical Loop Design
- **One kernel per (i,j) column, not per level:** Avoids redundant synchronization; allows stateful `tau_cum` accumulation within a single thread.
- **Grid bounds are crucial:** Always query `box.smallEnd(2)`, `box.bigEnd(2)`, never assume AMR resolution.

### Atmospheric State Access
- **Type-safe accessor:** Use `state(i, j, k, comp)` with known component indices (e.g., `Rho_comp`, `Temp_comp`).
- **Defensive clipping:** If temperature or density is unphysical (e.g., `T ≤ 0`), clip rather than crash; log to verbosity.

### GPU Memory & Reduction
- **Avoid device↔host copies in kernel:** Keep diagnostics on device until after the parallel loop completes.
- **Reduce-then-copy pattern:** Aggregate scalars in the kernel, copy once per time step.

### Documentation
- **Contracts first:** Define GPU-safety, grid-adaptivity, and I/O rules in markdown *before* coding.
- **Inline documentation:** Every function should have a docstring explaining what it does, expected input ranges, and any GPU-safety assumptions.

---

## References

- Toon et al., 1989: "Rapid calculation of radiative heating rates...", *J. Geophys. Res.*, 94, 16387–16405.
- Beer, A., 1852: "Bestimmung der Absorption des rothen Lichts in farbigen Flüssigkeiten", *Ann. Phys. Chem.*, 86, 78–88.
- AMReX GPU Guide: https://amrex-codes.github.io/amrex/docs_html/GPU.html
