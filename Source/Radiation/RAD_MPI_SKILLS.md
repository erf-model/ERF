# Radiation Module: MPI, GPU & Parallelization Skills

This document captures essential MPI, GPU, and AMReX parallelization patterns required for radiation module development. It serves as a reference for avoiding common pitfalls and understanding the architectural constraints that shape the radiation solver.

---

## Part A: GPU-Safe Kernel Design

### A.1 – Mark All Device Functions Properly

**Pattern:**
```cpp
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_flux(amrex::Real tau, amrex::Real T)
{
    // Device-safe code only
    return some_value;
}
```

**Why:**
- `AMREX_GPU_HOST_DEVICE` tells the compiler to generate both host and device versions
- `AMREX_FORCE_INLINE` improves performance by reducing function call overhead in kernels
- Without these, the function cannot be called from within a device lambda

**Common Mistake:**
```cpp
// ❌ WRONG: Missing GPU markers
amrex::Real compute_flux(amrex::Real tau, amrex::Real T) { ... }

// Later in kernel:
amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    f = compute_flux(tau, T);  // Compile error!
});
```

**Lesson:**
Every helper function intended for device-side use must have GPU markers from day one.

---

### A.2 – No Host-Side I/O in Device Lambdas

**Pattern:**
```cpp
// ✅ CORRECT: No I/O inside kernel
amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    amrex::Real flux = compute_flux(...);
    // Compute, don't print
});

// After kernel, on host:
if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Max flux: " << max_flux << "\n";
}
```

**Why:**
- `amrex::Print()` and file I/O are host-side only; calling them from device code is a runtime error
- On CPU/OpenMP, this might accidentally work (no device), masking the bug

**Common Mistake:**
```cpp
// ❌ WRONG: Printing inside kernel
amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    amrex::Real tau = ...;
    amrex::Print() << "tau = " << tau << "\n";  // Runtime error on GPU!
});
```

**Lesson:**
Reduce diagnostics (compute max/min/sum in kernel), copy back to host, print on host.

---

### A.3 – Reduce Before Copy: Device-Side Reduction Patterns

**Pattern:**
```cpp
amrex::Real max_tau_global = 0.0;

// On device: compute per-level max, reduce to global
amrex::Real max_tau_device = 0.0;
amrex::ParallelFor(bx, [=, &max_tau_device] AMREX_GPU_DEVICE (int i, int j, int k) {
    amrex::Real tau_local = ...;
    amrex::HostDevice::Atomic::Max(&max_tau_device, tau_local);
});

// Copy from device to host
amrex::Gpu::synchronize();
max_tau_global = max_tau_device;

// Print on host
if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Max tau: " << max_tau_global << "\n";
}
```

**Why:**
- Copying every data point is slow; copying one scalar is fast
- Device atomics allow in-kernel reduction without synchronization
- `amrex::Gpu::synchronize()` ensures all device work is complete before reading

**Lesson:**
Aggregate first, copy second, print third.

---

## Part B: Grid-Adaptive Vertical Integration

### B.1 – Query Grid Bounds from Box, Not Constants

**Pattern:**
```cpp
amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    int kmin = bx.smallEnd(2);
    int kmax = bx.bigEnd(2);
    
    for (int kk = kmin; kk <= kmax; ++kk) {
        // Vertical sweep
    }
});
```

**Why:**
- `bx.smallEnd(2)`, `bx.bigEnd(2)` are the actual grid bounds for this box
- Different AMR levels or domains may have different resolutions
- Hardcoded bounds cause silent bugs on coarse grids

**Common Mistake (Phase 1c Bug):**
```cpp
// ❌ WRONG: Hardcoded bounds
for (int k = 0; k < 50; ++k) {  // What if domain has only 32 levels?
    tau_cum += tau_per_layer;
}
```

**Lesson:**
Always derive loop bounds from the box or geometry; never hardcode them.

---

### B.2 – Vertical Sweep: One Thread per (i,j), Sequential k Loop

**Pattern:**
```cpp
// ✅ CORRECT: One kernel per (i,j) column
amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    // This runs once per (i,j); k is already iterated by AMReX
    // NO LOOP OVER k HERE!
});

// For a true vertical sweep (stateful tau_cum), use 2D iteration:
const auto& lo = bx.loVect();
const auto& hi = bx.hiVect();
amrex::ParallelFor(amrex::Box(lo[0], lo[1], 0, hi[0], hi[1], 0),
    [=] AMREX_GPU_DEVICE (int i, int j, int) {
        amrex::Real tau_cum = 0.0;
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            tau_cum += tau_per_layer;
            // Accumulate and compute fluxes
        }
    }
);
```

**Why:**
- Vertical integration is stateful (tau_cum accumulates); needs one thread per column
- Launching one kernel per level wastes thread blocks and requires thread-block-level reduction
- 2D parallelism (i, j) is natural; 1D k-loop is sequential within each thread

**Lesson:**
Design kernels around the problem structure: horizontal parallelism, vertical sequentiality.

---

### B.3 – Extract Box Bounds Correctly

**Pattern:**
```cpp
const auto& domain = geom[lev].Domain();
int kmin = domain.smallEnd(2);
int kmax = domain.bigEnd(2);

amrex::ParallelFor(amrex::Box(domain.smallEnd(0), domain.smallEnd(1), 0,
                              domain.bigEnd(0), domain.bigEnd(1), 0),
    [=] AMREX_GPU_DEVICE (int i, int j, int) {
        for (int k = kmin; k <= kmax; ++k) {
            // Vertical sweep
        }
    }
);
```

**Why:**
- `geom[lev].Domain()` is the full domain extent on this level
- `.smallEnd()`, `.bigEnd()` are 0-indexed inclusive bounds
- Creating a 2D `Box` with fixed `k` ensures parallel iteration only over (i, j)

**Lesson:**
Use geometry queries to make code robust to AMR refinement and domain size changes.

---

## Part C: Atmospheric State Access & Validation

### C.1 – Safe MultiFab Component Access

**Pattern:**
```cpp
// Define component indices (typically in ERF.H or a shared header)
constexpr int Rho_comp = 0;
constexpr int RhoU_comp = 1;
constexpr int RhoV_comp = 2;
constexpr int RhoW_comp = 3;
constexpr int RhoTheta_comp = 4;
constexpr int Temp_comp = 5;  // Diagnostic, computed from RhoTheta

// In kernel:
amrex::Real rho = state(i, j, k, Rho_comp);
amrex::Real T = state(i, j, k, Temp_comp);

// Use with defensive checks
if (rho <= 0.0) rho = 1.0;  // Clip to safe value
if (T <= 0.0) T = 288.15;   // Clip to sensible default
```

**Why:**
- Named component indices are more readable and less error-prone than magic numbers
- State variables come from prognostic arrays; defensive clipping prevents NaN/Inf
- Logging to verbosity helps diagnose physics errors

**Common Mistake:**
```cpp
// ❌ WRONG: Magic numbers, no validation
amrex::Real rho = state(i, j, k, 0);  // What is component 0?
amrex::Real T = rho / (287.0 * pressure);  // What if pressure is zero?
```

**Lesson:**
Use named constants, always validate ranges, log anomalies.

---

### C.2 – Defensive Clipping for Unphysical Values

**Pattern:**
```cpp
amrex::Real T = state(i, j, k, Temp_comp);
if (T <= 0.0 || T > 400.0) {
    T = 288.15;  // Clip to sensible default
    if (verbosity >= 1) {
        // Log once per simulation or to a counter
        // (avoid spamming if many points are bad)
    }
}

// Or for density:
amrex::Real rho = state(i, j, k, Rho_comp);
if (rho <= 0.0) {
    rho = 1.0;  // Typical sea-level density [kg/m^3]
}
```

**Why:**
- Simulation initialization or numerical schemes may produce unphysical intermediate values
- Crashing on the first bad value halts the run; clipping allows recovery
- Logging alerts developers to configuration issues

**Lesson:**
Be defensive in physics kernels; crash-on-error is better than silent NaN propagation, but clipping is better than both.

---

### C.3 – Logging Without Spamming

**Pattern:**
```cpp
if (verbosity >= 1 && amrex::ParallelDescriptor::IOProcessor()) {
    // Log once per time step, not per grid cell
    // Compute aggregate (max error, number of clipped points) on device
    // Copy to host, print once
    amrex::Print() << "Radiation: clipped " << n_clipped << " density values\n";
}
```

**Why:**
- Printing per grid cell generates millions of lines; unusable
- One message per time step or per level is informative
- `IOProcessor()` ensures only rank 0 prints (avoids duplicate output in MPI)

**Common Mistake:**
```cpp
// ❌ WRONG: Prints per grid cell
for (int i = ...; i <= ...; ++i) {
    for (int j = ...; j <= ...; ++j) {
        for (int k = ...; k <= ...; ++k) {
            if (T < 0) amrex::Print() << "Bad T at " << i << " " << j << " " << k << "\n";
        }
    }
}
// Output file becomes gigabytes!
```

**Lesson:**
Aggregate errors per level/time step; print summaries, not details.

---

## Part D: Known Issues & Workarounds

### D.1 – Phase 1c Bug: Hardcoded Grid Bounds (FIXED in Phase 2)

**Issue:**
Early Phase 1 prototypes used `for (int k = 0; k < 50; ++k)` instead of querying the actual domain. On a 32-level domain, this silently accessed out-of-bounds memory.

**Workaround (Phase 1c):**
Replaced hardcoded `50` with `geom[lev].Domain().length(2)`.

**Fix (Phase 2):**
All vertical loops now use `bx.smallEnd(2)`, `bx.bigEnd(2)` or domain queries.

**Lesson:**
Always derive loop bounds dynamically; hardcoded constants are a source of silent bugs.

---

### D.2 – Phase 1→2 Transition: Per-Column Kernel Design

**Issue:**
Attempting to launch one kernel per (i, j, k) level is tempting but inefficient: each level would need to synchronize, compute isolated flux values, and communicate back to host.

**Phase 2 Solution:**
Use 2D parallel iteration `(i, j)` with sequential k-loop inside lambda:
```cpp
const auto& domain = geom[lev].Domain();
amrex::ParallelFor(
    amrex::Box(domain.smallEnd(0), domain.smallEnd(1), 0,
               domain.bigEnd(0), domain.bigEnd(1), 0),
    [=] AMREX_GPU_DEVICE (int i, int j, int) {
        amrex::Real tau_cum = 0.0;
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            tau_cum += tau_per_layer;
            // Compute, reduce within one thread
        }
    }
);
```

**Lesson:**
Structure parallelism around the problem geometry: horizontal (embarrassingly parallel), vertical (inherently sequential).

---

### D.3 – Phase 2 Discovery: Temperature Conversion Complexity

**Issue:**
ERF stores `RhoTheta`, not `Theta`. Converting to temperature requires:
- Dividing by density: `Theta = RhoTheta / Rho`
- Accounting for Exner function: `T = Theta * (p / p_ref)^(R_d / cp)`

This requires access to pressure profile, which depends on model (hydrostatic, non-hydrostatic) and background state.

**Phase 2 Simplification:**
For Phase 2 diagnostics, use simplified conversion:
```cpp
amrex::Real theta = rho_theta / rho;
// Apply defensive clipping
theta = std::max(theta, 100.0);  // Minimum
theta = std::min(theta, 400.0);  // Maximum
```

**Future Fix (Phase 5+):**
When integrating heating into RhoTheta, use full Exner correction and background pressure profile.

**Lesson:**
Recognize when a simplification is acceptable (Phase 2: diagnostics only) vs. when full physics is needed (Phase 5+: coupling).

---

### D.4 – Temperature Not Directly Available on State

**Issue:**
ERF's state arrays typically store `RhoTheta`, not temperature. Computing `T` requires access to pressure, which depends on the model type and background state.

**Workaround:**
- Check if ERF provides a diagnostic `Temp_comp` component
- If not, compute locally: `T = (RhoTheta / Rho) * (p / p_ref)^(R_d / cp)`
- Cache the result to avoid recomputation per column

**Lesson:**
Know what's in the state arrays; don't assume standard atmospheric variable names.

---

### D.5 – Optical Depth Per Layer Can Vary (Future)

**Issue (Phase 3+):**
Phase 1/2 assume uniform `tau_per_layer`. Phase 3+ will vary optical depth with height (clouds, aerosols).

**Current Workaround:**
Use a simple parameter `tau_per_layer` for all levels; plan to generalize to `tau(k)`.

**Future Fix:**
Create a `tau_profile[k]` array, read per-level values in the vertical sweep.

**Lesson:**
Design loops to allow future per-level arrays; don't hardcode uniformity into the kernel structure.

---

### D.6 – Phase 2 New Issue: Limited State Access in Diagnostics Function

**Issue:**
The `compute_twostream_radiation_diagnostics()` function is called from `ERF_Advance.cpp` but doesn't have direct access to state MultiFab arrays.

**Current Workaround (Phase 2):**
Implemented the `vertical_two_stream_sweep()` kernel structure but cannot fully activate it without passing state as parameter.

**Future Fix (Phase 3+):**
Modify function signature to accept `const Vector<MultiFab>&` or pass a single-level state reference:
```cpp
void compute_twostream_radiation_diagnostics(
    int lev,
    int nstep,
    amrex::Real time_step,
    const amrex::MultiFab& state);  // Add this
```

**Lesson:**
When designing GPU kernels, ensure the calling context has access to required data structures before committing to the kernel design.

---

### D.7 – Phase 2b Discovery: Host Loop Calling Device Function (FIXED in Phase 2c)

**Issue:**
Phase 2b code had a GPU-safety bug: the device-side kernel `vertical_two_stream_sweep()` (marked `AMREX_GPU_DEVICE`) was called from a host-side nested loop:

```cpp
// ❌ WRONG: Host loop calling device function
for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
    for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        vertical_two_stream_sweep(i, j, bx, geom_lev, state_fab, rad_choice,
                                  max_heating_col, sw_flux_col, lw_net_col);
        // This is a runtime error on GPU; on CPU/OpenMP with tiling, it's a data race
        max_heating_global = std::max(max_heating_global, max_heating_col);
        sw_surface_sum += sw_flux_col;
    }
}
```

This violates GPU-safety rules in two ways:
1. Device functions can only be called from device code (e.g., inside a `ParallelFor` lambda)
2. Host-side accumulation (`+=`, `std::max()`) of results inside the loop causes data races when tiling is enabled

**Phase 2c Fix:**
Use `amrex::ParallelFor` with device-side reduction:

```cpp
// ✅ CORRECT: ParallelFor with device-side reduction
amrex::ParallelFor(xy_box,
    [=, &reduce_data] AMREX_GPU_DEVICE (int i, int j, int)
    {
        amrex::Real max_heating_col = 0.0;
        amrex::Real sw_flux_col = 0.0;
        amrex::Real lw_net_col = 0.0;
        
        vertical_two_stream_sweep(i, j, bx, geom_lev, state_arr, rad_choice,
                                  max_heating_col, sw_flux_col, lw_net_col);
        
        // Accumulate on device via ReduceOps, not on host
        reduce_data.join(max_heating_col, sw_flux_col, lw_net_col);
    }
);

// After kernel completes, copy reduced scalar from device to host
amrex::Gpu::synchronize();
auto reduce_tuple = reduce_data.value(reduce_ops);
max_heating_global = amrex::get<0>(reduce_tuple);
```

**Prevention Rule:**
A `grep` check to ensure kernels are called via `ParallelFor`:
```bash
grep -B3 "vertical_two_stream_sweep(" Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp | grep -q "ParallelFor" || echo "FAIL: kernel not called via ParallelFor"
```

**Lesson:**
Never call a device function from host-side code. Always use AMReX's parallel launch mechanisms (`ParallelFor`, `reduce`, etc.). Device functions can only be called from within device lambdas.

---

## Part E: Testing & Validation Checklist

### E.1 – GPU Compilation Check
- [ ] Code compiles with `CUDA_ARCH=all` (or target GPU arch)
- [ ] No warnings about uninitialized device memory
- [ ] No warnings about implicit host→device copies

### E.2 – Grid Adaptivity Check
- [ ] Run on coarse domain (e.g., 16×16×16 grid)
- [ ] Run on fine domain (e.g., 256×256×256 grid)
- [ ] Verify vertical loop counts match expected resolution
- [ ] Check CSV output reflects correct layer count

### E.3 – Physical Correctness Check
- [ ] SW flux decreases monotonically from TOA to surface (for increasing optical depth)
- [ ] LW flux changes smoothly with temperature profile
- [ ] Heating rates have correct sign and reasonable magnitude
- [ ] No NaN/Inf in output (run with `IEEE=1` to catch NaN earlier)

### E.4 – State Variable Check
- [ ] Temperature and density read correctly
- [ ] No clipping errors in typical runs (verbosity log clean)
- [ ] Defensive clipping activates only in extreme cases

### E.5 – Performance Check
- [ ] Kernel time doesn't scale poorly with grid size
- [ ] No excessive host↔device transfers
- [ ] Reduction overhead is negligible (~1% of total time)

---

## Part F: Documentation Standards

### F.1 – Function Docstrings

**Pattern:**
```cpp
/**
 * @brief Compute direct-beam solar flux using Beer-Lambert law.
 *
 * Implements the formula:
 *   F_dir(z) = S0 * cos(zenith) * exp(-tau_cumulative(z) / cos(zenith))
 *
 * This is a GPU-safe inline function intended for use in device-side kernels.
 *
 * @param[in] tau_cumulative Cumulative optical depth from TOA [unitless].
 * @param[in] S0 Solar constant [W/m^2]. Must be positive.
 * @param[in] cos_zenith Cosine of solar zenith angle [unitless, 0 ≤ cos_zenith ≤ 1].
 *
 * @return Direct-beam flux [W/m^2]. Non-negative.
 *
 * @note If cos_zenith ≤ 0 (night), returns 0.
 * @note This function does NOT handle scattering or diffuse radiation (Phase 1 simplification).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_sw_direct_flux(amrex::Real tau_cumulative, amrex::Real S0,
                                    amrex::Real cos_zenith);
```

### F.2 – File Header

**Pattern:**
```cpp
/**
 * @file ERF_TwoStreamSW.H
 * @brief Clear-sky shortwave radiation using Beer-Lambert direct-beam model.
 *
 * Implements a simplified, clear-sky shortwave radiation model
 * using the Beer-Lambert law. Supports both Phase 1 (uniform optical depth)
 * and Phase 3+ (vertical optical depth profiles).
 *
 * All functions are GPU-safe and intended for use in device-side kernels.
 *
 * References:
 * - Beer, A., 1852: ...
 * - Bird et al., 1984: ...
 */
```

### F.3 – Code Comments

**Pattern:**
```cpp
// Vertical sweep: accumulate optical depth from TOA downward
amrex::Real tau_cum = 0.0;
for (int k = kmin; k <= kmax; ++k) {
    tau_cum += tau_per_layer;  // Assume uniform for Phase 1/2
    
    // Compute direct-beam flux at this level
    amrex::Real F_dir = compute_sw_direct_flux(tau_cum, S0, cos_zenith);
    
    // Store or accumulate heating (Phase 2: diagnostics only)
    diagnostics[i][j][k] = F_dir;
}
```

**Why:**
- Comments explain *why*, not *what* (code shows what)
- Phase labels (Phase 1/2, etc.) guide future maintainers
- References to formulas help with code review

---

## Summary: The Five Sacred Rules

1. **Mark functions:** `AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE` for device use
2. **No host I/O in kernels:** Reduce, copy, print in that order
3. **Dynamic bounds:** Query box/domain, never hardcode loop limits
4. **Per-(i,j) kernels:** Vertical sweep is sequential, horizontal is parallel
5. **Defensive clipping:** Validate state, clip unphysical values, log aggregates

---

## References

- AMReX GPU Guide: https://amrex-codes.github.io/amrex/docs_html/GPU.html
- AMReX Parallel Loop Patterns: https://amrex-codes.github.io/amrex/docs_html/GPU_HowTo.html#parallel-for
- Atomic Operations: https://amrex-codes.github.io/amrex/docs_html/GPU_HowTo.html#atomic-operations
