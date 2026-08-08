# Radiation Module: MPI, GPU & Parallelization Skills

This document captures essential MPI, GPU, and AMReX parallelization patterns required for radiation module development. It serves as a reference for avoiding common pitfalls and understanding the design rationale behind the two-stream radiation implementation.

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

**Status Update (Phase 3):** Implemented via `tau_layer_value()`, which
adds a cloud-layer contribution based on height rather than a fully
general `tau(k)` array — sufficient for the height-band cloud model but
still not a fully arbitrary per-level profile.

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

### D.8 – Phase 3 Lesson: Re-Verify Every `_enabled`-Style Flag After Any Rewrite

**Issue:**
Phase 2d discovered that `sw_enabled` / `lw_enabled` gating had silently stopped being respected in some code paths after an earlier rewrite of `vertical_two_stream_sweep()` and `compute_twostream_radiation_diagnostics()`.

**Prevention Rule:**
After ANY rewrite of a function containing `_enabled`-style (or enum-based, e.g., `tau_profile_type`, `isothermal_test`) gating flags, explicitly grep the modified file for every place the flag SHOULD gate behavior, and confirm each one still does:

```bash
grep -n "sw_enabled\|lw_enabled\|isothermal_test\|tau_profile_type\|cloud_fraction" \
  Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp
```

Then manually confirm, for each match, that it actually wraps the corresponding SW/LW computation rather than merely appearing in an unrelated comment or unused branch.

**Lesson:**
A boolean flag with no gating effect is a silent physics bug, not a compile error. Every `_enabled`/`_type`-style ParmParse flag must be re-verified — not just re-read — after any rewrite of the function(s) it's supposed to gate.

---

### D.9 – Phase 3 Lesson: A Task Is Not Complete Until It Is Pushed and a PR Exists

**Issue:**
Prior attempts at delivering radiation-module phases have produced fully-correct, well-documented code changes that were never actually pushed to a remote branch, or where a PR was never opened despite the branch existing.

**Prevention Rule:**
- Never use `git branch -r` alone to conclude a remote branch does not exist — this only lists branches your local repo already has a tracking ref for. Use `git ls-remote origin <branch>` or an equivalent server-side check instead.
- After pushing, verify success by checking actual command output for a new/updated ref (not just a zero exit code), and cross-check with `git ls-remote origin | grep <branch-name>`.
- If a push fails, retry exactly once. If it fails again, stop and report the exact error message verbatim — do not attempt further workarounds, and do not describe the task as complete.
- A local commit with no pushed branch, or a pushed branch with no opened pull request, is an **incomplete task**, regardless of how complete or well-documented the code changes themselves are.

**Lesson:**
Code correctness and delivery mechanics (push + PR) are two independent completion criteria. Verify both explicitly before ending a session.

---

### D.10 – Phase 4 Lesson: A Stale Hardcoded Phase Tag in Shared Diagnostics Code Is a Silent, Cross-Phase Bug

**Issue:**
`RadiationDiagnostics::append()` (`Source/Radiation/ERF_RadiationDiagnostics.cpp`) hardcoded its bracketed debug tag as `[RAD][Phase1][RadiationDiagnostics::append]` from Phase 1 onward. Because this diagnostics class is shared and reused unchanged across every subsequent phase (2, 2b, 2c, 2d, 3, 4, ...), the tag was never updated — every debug print, in every phase, silently claimed to be "Phase1" output. This was discovered when running the Phase 4 `SW_Scattering_Cloud` RegTest with `erf.radiation.v = 1`: the log clearly showed Phase 4 cloud/scattering diagnostics under a `[Phase1]` tag, which would mislead anyone grepping logs by phase to debug a specific phase's behavior.

**Root Cause:**
Any hardcoded phase/version label embedded in a **shared, cross-phase utility class** (as opposed to phase-specific driver code like `vertical_two_stream_sweep()`) will silently go stale the moment a new phase starts using that utility, because there is no compiler or test failure to catch it — it's a purely cosmetic/logging bug with no functional effect on flux correctness.

**Fix:**
Removed the hardcoded `[Phase1]` segment; the tag is now the phase-agnostic `[RAD][RadiationDiagnostics::append]`, which remains accurate regardless of which phase's code is calling it. The always-on `RADIATION_DIAG:` line used by RegTest check scripts (`check_flux_accuracy.py` etc.) was left untouched, since it never included a phase tag and is unaffected.

**Prevention Rule:**
When adding a new phase, grep all touched (and reused-but-unmodified) files for hardcoded phase-number strings, not just files being actively edited:
```bash
grep -rn '\[Phase[0-9]\]' Source/Radiation/
```
Any hit inside a file that is *shared* across phases (i.e., not itself part of the phase-specific driver logic being changed) should be treated as a latent staleness bug, even if that file isn't otherwise part of the current phase's diff.

**Lesson:**
A hardcoded phase label is fine in phase-specific narrative comments/documentation (where it's historically accurate and never changes), but must NEVER appear in a runtime string embedded in a class/function that is reused, unmodified, across multiple phases. Prefer phase-agnostic tags (`[RAD][...]`) for anything printed by shared infrastructure code.

---

### D.11 – Phase 5 Lesson: "Correct Code" and "Wired-In Code" Are Different Completion Criteria

**Issue:**
Phases 1-4 of the two-stream radiation module were fully implemented,
documented, hand-traced, and RegTested — but a repo-wide search at the
start of Phase 5 confirmed `compute_twostream_radiation_diagnostics()`
was never called from anywhere in the codebase. `advance_radiation()`
(`Source/TimeIntegration/ERF_AdvanceRadiation.cpp`) only drove the
RRTMGP path (`solverChoice.rad_type != RadiationType::None`), which uses
a completely different enum from the TwoStream path
(`solverChoice.radChoice.rad_type == RadType::TwoStream`). Additionally,
`qheating_rates[lev]` was only allocated for the RRTMGP path, so even if
the TwoStream driver had been called, it would have had nowhere to write
its heating rates, and the RhoTheta source-term injection in
`ERF_MakeSources.cpp` was similarly gated on RRTMGP only.

**Root Cause:**
Four independent, all-or-nothing gates (allocation, driver call, source
injection, and the underlying per-level output the driver produces) all
had to individually recognize the new radiation-model enum value. Adding
correct physics to a driver function does not automatically make that
function get called, and even if called, its output MultiFab might not
exist, and even if it exists and is populated, the downstream consumer
might not read from it for this code path.

**Prevention Rule:**
Before considering any new physics module "complete," explicitly verify —
via code search, not just code review of the files you touched — that:
1. The new driver function is actually invoked somewhere in the
   simulation's time loop (`grep -rn "function_name(" Source/` and
   manually confirm at least one call site is unconditionally reachable
   under the intended `ParmParse` configuration).
2. Any output data structure it writes into is allocated under the same
   condition that gates the call (not a different, historically-earlier
   condition that happened to cover only the previous physics option).
3. Any downstream consumer of that output data structure is gated to
   include the new physics option, not just the original one(s) it was
   written for.

**Lesson:**
"The math is right and the RegTest for the isolated function passes" and
"this function affects the simulation" are two independent claims. This
codebase-wide, four-gate discovery (allocation → driver call → per-level
computation → source injection) — where physics work from four prior
phases turned out to be entirely inert — is the strongest possible
illustration of why the D.9 delivery-mechanics lesson (push + PR are
independent of code correctness) generalizes: *wiring* is also
independent of *correctness*, and must be verified with the same rigor.

---

### D.12 – Phase 6 Lesson: Temporal Wiring Errors Hide as Duplicate Rows or Missing Forcing

**Issue:**
Phase 5 wired `advance_radiation()` into the time loop, successfully computing
and injecting heating rates. However, without explicit auditing of *when*
the function is called relative to the slow/fast substep structure, a subtle
temporal error could occur:

- If `advance_radiation()` were called multiple times per slow step (e.g.,
  once before dycore and once after), or if the `is_slow_step` gate in
  `ERF_MakeSources.cpp` were incorrect (e.g., applied on fast steps when
  it should only apply to slow steps), the radiation forcing could be applied
  multiple times per physical slow step, resulting in:
  - Duplicate diagnostic CSV rows per timestep (confusing row counts)
  - Unintended effective forcing multiplication (2x, 3x the intended heating)
  - Numerical instability if heating rates are large

- Conversely, if the gate is too restrictive (e.g., gated on fast-step-only),
  no radiation forcing would occur at all, and the heating-rate diagnostics
  might still populate the CSV (making the wiring appear active) even though
  the source term never actually affects the prognostic state.

**Root Cause:**
A single driver function call (`compute_twostream_radiation_diagnostics()`)
paired with a single source-term gate (`is_slow_step`) creates TWO independent
control points that must remain synchronized. Auditing that they are actually
in sync requires tracing through:
1. Where and how often `advance_radiation()` is called in the time loop
2. Where and when `is_slow_step` becomes true/false during the step sequence
3. Whether diagnostic CSV logging happens once or multiple times per step

This is not primarily a *coding* error (like a missing GPU marker) but a
*temporal/architectural* error that manifests as unexpected behavior in output.

**Prevention Rule:**
When wiring a new physics driver into the time loop (especially one with
both diagnostic output and source-term injection), document and audit:

1. **Call Site Audit:** Grep for where the driver is called and verify
   manually that it is called exactly as many times per slow step as intended:
   ```bash
   grep -rn "advance_radiation\|compute_twostream_radiation" Source/
   ```
   Then trace each call site to confirm it is in the intended loop/condition.

2. **Source-Term Gate Audit:** Find where source terms from this driver are
   injected and verify the gate matches the call site's temporal position:
   ```bash
   grep -rn "is_slow_step.*qheating\|qheating.*is_slow_step" Source/
   ```
   Confirm that if the driver is called once per slow step, the gate is
   `is_slow_step==true` (and vice versa).

3. **Diagnostic Output Cadence:** Check the diagnostic CSV file format and
   expected row count:
   - If running for N slow steps, expect approximately N rows (or N × M if
     the driver is called M times per slow step by design).
   - If the row count is 2x expected or differs unexpectedly, suspect either:
     - The driver is being called twice per step
     - The diagnostic logging is happening in two places
     - The source-term gate is asymmetric relative to the driver call

4. **Cross-Step Consistency:** For multi-step runs, verify that diagnostic
   values do not have unintended jumps or discontinuities that would suggest
   the source term is being applied inconsistently across steps.

**Phase 6 Audit Findings (Example):**
For the TwoStream two-stream radiation:
- `advance_radiation()` is called once per `ERF::Advance()`, which is once
  per slow step (line 150 in `ERF_Advance.cpp`, before dycore)
- `qheating_rates[lev]` is gated on `is_slow_step` in `ERF_MakeSources.cpp`
  (line 268)
- Diagnostic CSV output is appended once per call to
  `compute_twostream_radiation_diagnostics()`
- Expected row count is 1 row per slow step (verified via Phase 5/6 RegTests)

**Lesson:**
Temporal wiring errors do not produce compile errors or obvious NaN values.
Instead, they produce subtle inconsistencies: extra diagnostic rows, apparent
failure of source terms to affect the solution, or doubling of effective
forcing. Always audit temporal semantics explicitly; do not assume single
calls or correct gating without tracing the code paths involved.

### D.13 – Phase 6 Lesson: Treat Diagnostics Identity as `(step, time, call_site)`, Not `step` Alone

**Issue:**  
Phase 6 introduced intentional multi-call-site diagnostics per coarse step
(`pre_dycore`, `post_dycore`). A legacy duplicate guard keyed only on `step`
silently dropped valid rows, while checker logic keyed only on `(step,time)`
misinterpreted valid entries as duplicates.

**Root Cause:**  
Diagnostics cadence evolved (2 rows/step by design), but identity semantics in
both writer and checker were still based on earlier single-row assumptions.

**Fix:**
1. Extend diagnostics schema with `call_site`.
2. Pass `call_site` through `compute_twostream_radiation_diagnostics(...)`
   and into `RadiationDiagnostics::append(...)`.
3. Update duplicate guard to preserve valid pre/post rows while rejecting only
   true duplicates of the same event.
4. Update checker expectations to:
   - expect 2 rows/step in Phase 6 timing test
   - validate one `pre_dycore` and one `post_dycore` per step
   - treat `(step,time,call_site)` as uniqueness key

**Recommended guard pattern:**
```cpp
if (step == m_last_write_step &&
    call_site == m_last_write_call_site &&
    std::abs(time - m_last_write_time) < 1.0e-12) {
  return;
}
m_last_write_step = step;
m_last_write_call_site = call_site;
m_last_write_time = time;
```

**Prevention Rule:**  
Any time diagnostics cadence changes (new phase hooks, pre/post hooks, nested
driver calls), update all three together:
1. Writer schema/fields
2. Duplicate guard identity key
3. Checker row-count and uniqueness logic

**Lesson:**  
A cadence change without an identity-model update creates “false duplicate”
bugs that look like physics regressions but are actually observability defects.

---

### D.14 – Phase 6 Lesson: Don’t Infer Call Site from Out-of-Scope Time Variables

**Issue:**  
An attempted `call_site` inference in `ERF_AdvanceTwoStreamRadiation.cpp` used undeclared variables (e.g., `old_time`, `dt_lev`) and failed compilation.

**Root Cause:**  
Call-site classification was attempted inside a function that did not own the full time-loop context.

**Fix:**  
Thread `call_site` explicitly from the caller (time-integration layer) into:
- `compute_twostream_radiation_diagnostics(...)`
- `RadiationDiagnostics::append(...)`

**Prevention Rule:**  
If a function lacks authoritative temporal context, pass semantic context as an explicit argument rather than reconstructing it indirectly.

**Lesson:**  
Semantic labels (`pre_dycore` / `post_dycore`) belong to call-site ownership, not heuristic inference inside lower-level diagnostics code.

---

### D.15 – Phase 6 Lesson: Checker File Path Must Follow Runtime `diag_file`

**Issue:**  
Phase 6 checker initially looked for `radiation_diagnostics.csv`, while the run wrote `radiation_phase6_timing_diag.dat`, causing a false failure (`file not found`).

**Root Cause:**  
Checker used a hardcoded filename instead of following configured output path.

**Fix:**  
Read the configured diagnostics filename (or support expected candidates), and document the active file in test output.

**Prevention Rule:**  
Never hardcode diagnostics filename in validation scripts when runtime config controls output naming.

**Lesson:**  
Path mismatches can masquerade as model failures; validate I/O contract first.

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
- [ ] (Phase 4+) Diffuse SW flux is exactly zero when `single_scattering_albedo`/`cloud_single_scattering_albedo` are 0.0, and strictly positive when nonzero and a direct beam is present
- [ ] (Phase 5+) New driver functions are confirmed CALLED from the time
      loop via code search (not just present in the source tree), their
      output MultiFabs confirmed ALLOCATED under the matching gate, and
      downstream consumers confirmed to READ from them under that same
      gate — see D.11

### E.4 – State Variable Check
- [ ] Temperature and density read correctly
- [ ] No clipping errors in typical runs (verbosity log clean)
- [ ] Defensive clipping activates only in extreme cases

### E.5 – Performance Check
- [ ] Kernel time doesn't scale poorly with grid size
- [ ] No excessive host↔device transfers
- [ ] Reduction overhead is negligible (~1% of total time)

### E.6 – Cross-Phase Regression Check (Phase 4+)
- [ ] Grep all shared/reused files (not just newly-modified ones) for
  hardcoded phase-number labels (`grep -rn '\[Phase[0-9]\]'
  Source/Radiation/`) — see D.10
- [ ] Re-run all prior phases' RegTests (not just the new phase's
  RegTest) and confirm byte-identical or numerically-identical output
  to the previous phase, per each new phase's own hand-traced arithmetic

### E.7 – Phase 6 Diagnostics Cadence Validation (add to checklist)

- [ ] Diagnostics schema includes `call_site` column
- [ ] `RADIATION_DIAG:` log lines include `call_site`
- [ ] Phase 6 timing test expects 2 rows per coarse step
- [ ] Each step has one `pre_dycore` and one `post_dycore`
- [ ] Duplicate guard preserves distinct pre/post rows
- [ ] Checker reads configured diagnostics file path (not hardcoded default)

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
 * @note This function does NOT handle scattering or diffuse radiation (Phase 1 simplification;
 *       see compute_sw_diffuse_flux() added in Phase 4 for the diffuse component).
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
 * @brief Clear-sky shortwave radiation using Beer-Lambert direct-beam model,
 * plus (Phase 4) a two-stream diffuse (scattering) approximation.
 *
 * Implements a simplified, clear-sky shortwave radiation model
 * using the Beer-Lambert law. Supports both Phase 1 (uniform optical depth)
 * and Phase 3+ (vertical optical depth profiles). Phase 4 adds a diffuse
 * SW flux term via the Meador-Weaver (1980) two-stream approximation.
 *
 * All functions are GPU-safe and intended for use in device-side kernels.
 *
 * References:
 * - Beer, A., 1852: ...
 * - Bird et al., 1984: ...
 * - Meador and Weaver, 1980: ...
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

**Note (Phase 4 addendum, see D.10):** Phase labels in *comments and
documentation* (like the ones above) are fine and encouraged — they are
historically accurate and never need to change. The bug D.10 describes is
specifically about hardcoded phase labels embedded in **runtime debug
print strings** within shared, reused-across-phases code, which silently
goes stale. Keep these two categories distinct.

---

## Summary: The Eight Sacred Rules (append)

1. **Mark functions:** `AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE` for device use  
2. **No host I/O in kernels:** Reduce, copy, print in that order  
3. **Dynamic bounds:** Query box/domain, never hardcode loop limits  
4. **Per-(i,j) kernels:** Vertical sweep is sequential, horizontal is parallel  
5. **Defensive clipping:** Validate state, clip unphysical values, log aggregates  
6. **No stale phase tags in shared runtime strings:** Avoid hardcoded `[Phase#]` in cross-phase utility logs  
7. **Wiring is independent of correctness:** Verify call path, allocation gate, and downstream consumption all match  
8. **Phase 6 cadence rule:** When diagnostics cadence includes multiple call-sites per step, treat event identity as `(step, time, call_site)`, pass `call_site` explicitly from time-loop call sites, and keep checker file path/row-count logic aligned with runtime configuration.

---

## D.16 – Lesson (Phase 7): Diagnostics Controls Must Not Alter Physics Path

**Goal:**
Phase 7 introduces runtime controls for diagnostics output cadence and streams.
This lesson ensures diagnostics controls remain orthogonal to the physics
solver and do not leak into heating calculations or source-term application.

**Common Mistake:**
```cpp
// ❌ WRONG: Letting diag_callsite_mode affect physics
void compute_twostream_radiation(...) {
    // Physics computation
    if (diag_callsite_mode == "pre_only") {
        // Only compute heating during pre-step
        // WRONG: This changes the physics!
    }
}
```

**Correct Pattern:**
```cpp
// ✅ CORRECT: Diagnostics gating only affects output, not computation
void compute_twostream_radiation(...) {
    // Compute heating unconditionally
    compute_sw_heating(qheating_rates, ...);
    compute_lw_heating(qheating_rates, ...);
    
    // Only the diagnostics append respects callsite mode:
    RadiationDiagnostics rad_diag(...);  // controls passed here
    rad_diag.append(...);  // filters output only
}
```

**Why:**
- Diagnostics controls are for observability; physics must be deterministic
- Heating rates go into source terms regardless of output mode
- Reproducibility requires decoupled layers: compute, then observe
- Users configuring diagnostics must not inadvertently change results

**Phase 7 Validation Checklist:**
- [ ] All `diag_*` parameters used **only** in `RadiationDiagnostics` append path
- [ ] No `diag_*` check in physics solver (`vertical_two_stream_sweep`, heating rate computation)
- [ ] No `diag_*` check in source-term injection code
- [ ] Heating rates computed and stored identically regardless of `diag_callsite_mode`
- [ ] CSV/stdout enable flags gate **output only**, not computation

---

## D.17 – Lesson (Phase 8): Benchmark Reproducibility & Diagnostics-Mode-Aware Validation

**Goal:**
Phase 8 establishes a canonical benchmark suite for TwoStream radiation. This lesson
ensures benchmarks are reproducible and validates that Phase 6/7 diagnostics controls
work correctly without affecting physics.

**Key Challenge:**
Validation metrics (flux, heating, cadence) are extracted from diagnostic CSV output.
The CSV structure depends on runtime configuration (especially `diag_callsite_mode`).
A robust validation suite must:
1. Be agnostic to physics changes (validate existing behavior, not drive new physics)
2. Account for diagnostics mode when checking row counts and call-site presence
3. Provide deterministic pass/fail criteria via centralized tolerance config
4. Support CI/CD automation with clear exit codes (0=pass, 1=fail)

**Common Mistake:**
```cpp
// ❌ WRONG: Benchmark metrics hard-coded to specific mode
def check_row_count(csv_path, expected_rows):
    rows = read_csv(csv_path)
    if len(rows) != expected_rows:
        fail()  # Fails if mode changes from "both" to "pre_only"!
```

**Correct Pattern:**
```cpp
// ✅ CORRECT: Metrics aware of configured diagnostics mode
class BenchmarkCase:
    diag_callsite_mode: str          # "both", "pre_only", or "post_only"
    rows_per_step: int               # Derived from mode
    expected_diag_rows: int          # = rows_per_step * num_steps

def validate_case(case, csv_path):
    rows = read_csv(csv_path)
    expected = case.expected_diag_rows
    
    # Check total row count (±tolerance for startup/teardown)
    if abs(len(rows) - expected) > ROW_COUNT_ABS_TOL:
        fail(f"Expected {expected}±{TOL}, got {len(rows)}")
    
    # Check per-step multiplicity
    step_counts = Counter(step for row in rows)
    if any(count != case.rows_per_step for count in step_counts.values()):
        fail("Step multiplicity doesn't match mode")
    
    # Validate call-site filtering
    if case.diag_callsite_mode == "pre_only":
        if not all("pre" in row["call_site"] for row in rows):
            fail("Expected only pre call-sites for pre_only mode")
    elif case.diag_callsite_mode == "both":
        if not (any("pre" in row["call_site"] for row in rows) and
                any("post" in row["call_site"] for row in rows)):
            fail("Expected both pre and post call-sites for both mode")
```

**Why:**
- Benchmark cases define their diagnostics configuration (mode, frequency)
- CSV schema depends on this configuration (row count, call_site tags)
- Tolerances must be physics-independent (pure metric checks)
- Validation suite must be robust to configuration changes
- Decoupling configuration (BenchmarkCase) from checks (validate_case()) keeps 
  code maintainable as diagnostics controls evolve

**Phase 8 Validation Checklist:**
- [ ] Each benchmark case defines: `diag_callsite_mode`, `rows_per_step`, `expected_diag_rows`
- [ ] Row count checks use case-specific expected values (not hard-coded)
- [ ] Call-site validation respects mode ("both" checks for both, "pre_only" checks for pre only)
- [ ] Metrics extracted are **independent of physics** (purely from CSV)
- [ ] Tolerances centralized in `benchmark_tolerances.py` (no magic numbers in checks)
- [ ] Pass/fail logic deterministic and CI-friendly (non-zero exit code on any fail)
- [ ] Reports machine-readable (JSON) and human-readable (Markdown)
- [ ] Phase 6/7 "both" mode case passes with Phase 7 `check_timing_consistency.py`
- [ ] Phase 8 "pre_only" mode case passes with mode-aware row count and call-site checks

**Phase 8 Reproducibility Requirements:**
- Fixed `dt` and `stop_time` per case (deterministic step count)
- Explicit `diag_file` path per case (no default inference)
- Explicit `diag_callsite_mode` per case (documented in case definition)
- Deterministic case execution order (sorted by case name)
- Deterministic report row order (sorted by case name)
- Seed/RNG settings (if applicable) documented
- Metric tolerances stored in version-controlled config

**Integration with CI/CD:**
```bash
#!/bin/bash
# Run all 5 benchmark cases
cd Exec/CanonicalTests/Radiation
for case in SW_ClearSky_Analytical LW_Isothermal SW_Cloud_Layer SW_Scattering_Cloud Phase8_Benchmark_Suite/cases/phase6_timing; do
    cd $case && mpirun -np 1 erf.ex inputs || exit 1
    cd -
done

# Validate all outputs
cd Phase8_Benchmark_Suite
python3 run_benchmark_suite.py --verbose
# Exit code 0 = all pass; 1 = any fail
exit $?
```

---

## References

- AMReX GPU Guide: https://amrex-codes.github.io/amrex/docs_html/GPU.html
- AMReX Parallel Loop Patterns: https://amrex-codes.github.io/amrex/docs_html/GPU_HowTo.html#parallel-for
- Atomic Operations: https://amrex-codes.github.io/amrex/docs_html/GPU_HowTo.html#atomic-operations
- Meador, W. E., and W. R. Weaver, 1980: "Two-stream approximations to radiative transfer in planetary atmospheres", *J. Atmos. Sci.*, 37, 630–643.

---

## Phase 9: Diagnostics Dedup & Nonuniform dz Lessons Learned

### Lesson 1: Triple-Tuple Dedup Identity is Robust

**Problem**: Early diagnostics logic used only (step, time) for dedup, which accidentally collapsed legitimate multi-call-site entries (e.g., pre_dycore and post_dycore at the same step).

**Solution (Phase 9)**: Use strict 3-tuple identity:
```
(step, call_site, time)
```

This ensures:
- Accidental repeated calls to append() at identical (step, call_site, time) are suppressed.
- Legitimate pre + post entries differ in call_site, so identity tuple differs → both retained.
- Mode filtering (pre_only, post_only, both) is orthogonal to dedup.

**Key Implementation Detail**:
```cpp
// Phase 9: Dedup guard order is critical
// 1. First: check mode filtering (returns early if unwanted)
// 2. Then: check (step, call_site, time) dedup tuple
// 3. Update: save last_step, last_call_site, last_time for next call

if (!should_emit_for_this_mode) return;  // Mode filter first
if (step == m_last_write_step && call_site == m_last_write_call_site && 
    time close to m_last_write_time) return;  // Then dedup
// Update guard only if we're about to write
m_last_write_step = step;
m_last_write_call_site = call_site;
m_last_write_time = time;
```

**Defensive Implication**:
- Time tolerance (m_diag_dedup_tol) is needed because multiple append() calls at slightly different times (FP rounding) should be treated as "same" event.
- Call_site must be exact string match; partial matching (e.g., regex) can mask bugs.

### Lesson 2: Nonuniform dz Framework Must Have Uniform Fallback

**Problem**: Heating divergence formulas divide by dz, which must be physically accurate for terrain-aware grids, but uniform grids benefit from algorithmic simplicity and performance.

**Solution (Phase 9)**: 
- Local per-level dz array `dz_level[MAX_RAD_LEVELS]`
- Initialize all entries to uniform `geom.CellSize(2)` (current behavior)
- When terrain/nonuniform support is added, populate from z_phys_cc or equivalent
- Defensive fallback: bounds-check array access; return uniform dz if out-of-bounds

```cpp
// Phase 9: Per-level dz framework
amrex::Real dz_level[MAX_RAD_LEVELS];
for (int k = 0; k < nlev; ++k) {
    dz_level[k] = dz_uniform;  // Current: always uniform
    // Future: dz_level[k] = z_cc[k] - z_cc[k+1]  (when z_phys_cc available)
}

// At heating calculation:
int k_idx = k - kmin;
amrex::Real dz_heating = (k_idx >= 0 && k_idx < MAX_RAD_LEVELS) 
    ? dz_level[k_idx] 
    : dz_uniform;  // Fallback
amrex::Real Q_sw = compute_sw_heating_rate(..., dz_heating, ...);
```

**Key Insight**:
- Separating uniform dz (for cloud-layer height) from per-level dz (for heating) keeps cloud detection logic unchanged.
- Cloud layers are currently defined in height space (m), not index space, so they should use uniform dz for now.
- Future: When cloud detection becomes terrain-aware, use z_phys_cc height directly (avoid dz altogether).

**Defensive Implication**:
- Heating rate functions must guard against dz <= 0 (returns 0, not NaN).
- Edge cases: Very thin layers (dz ~1 cm) or thick layers (dz ~1 km) should both work correctly.
- Extreme flux divergence (e.g., tau very large) can produce large heating rates; guard with isfinite().

### Lesson 3: Finite Checks Prevent Silent Corruption

**Problem**: If dz becomes zero or negative (grid bug), the heating rate is NaN. If qheating_rates gets populated with NaN, it silently corrupts RhoTheta in later source-term injection.

**Solution (Phase 9)**:
1. Input guards in `compute_sw_heating_rate()` / `compute_lw_heating_rate()`:
   ```cpp
   if (dz <= 0.0 || rho <= 0.0 || cp <= 0.0) return 0.0;
   ```

2. Output guards:
   ```cpp
   amrex::Real heating = flux_divergence / (rho * cp);
   if (!std::isfinite(heating)) return 0.0;  // Catch NaN/Inf
   ```

3. No warning/error log: Silently return 0 to avoid log spam in large simulations. (Future: Phase 10 can add counters to track how often this happens.)

**Defensive Implication**:
- These guards are GPU-safe inline functions.
- Returning 0 heating rate is physically conservative (assumes no radiative forcing, which is safe if inputs are invalid).
- Silent suppression is acceptable for an integration-polish phase; explicit warnings/assertions can be added later if edge cases become common.

### Phase 9 Takeaway

**Dedup + Nonuniform dz + Finite checks** = robust integration polish that:
- Prevents silent data corruption (dedup tuples, finite checks)
- Prepares for terrain-aware grids without breaking uniform-grid behavior
- Remains GPU-safe and performance-neutral on current (uniform) simulations
