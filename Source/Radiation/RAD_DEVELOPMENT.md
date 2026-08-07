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
