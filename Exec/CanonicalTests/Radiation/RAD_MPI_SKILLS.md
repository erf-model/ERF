# Radiation Module: MPI, GPU, and Parallelization Skills

This document captures implementation practices for developing and maintaining the radiation module within ERF's MPI- and GPU-enabled execution model. It preserves the original technical lessons while reorganizing them into stable engineering guidance for kernel design, data access, diagnostics, validation, and documentation.

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
- On CPU/OpenMP, this might accidentally work, masking the bug

**Common Mistake:**
```cpp
// ❌ WRONG: Printing inside kernel
amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    amrex::Real tau = ...;
    amrex::Print() << "tau = " << tau << "\n";
});
```

**Lesson:**
Reduce diagnostics on device, transfer only reduced values, and print on host.

---

### A.3 – Reduce Before Copy: Device-Side Reduction Patterns

**Pattern:**
```cpp
amrex::Real max_tau_global = 0.0;

amrex::Real max_tau_device = 0.0;
amrex::ParallelFor(bx, [=, &max_tau_device] AMREX_GPU_DEVICE (int i, int j, int k) {
    amrex::Real tau_local = ...;
    amrex::HostDevice::Atomic::Max(&max_tau_device, tau_local);
});

amrex::Gpu::synchronize();
max_tau_global = max_tau_device;

if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Max tau: " << max_tau_global << "\n";
}
```

**Why:**
- Copying every data point is slow; copying one scalar is fast
- Device atomics allow in-kernel reduction without synchronization
- `amrex::Gpu::synchronize()` ensures device work is complete before reading

**Lesson:**
Aggregate first, copy second, print third.

---

### A.4 – Shared Runtime Strings Must Be Phase-Agnostic

**Issue:**
A shared diagnostics utility once hardcoded `[Phase1]` in runtime debug output. Because the same utility was reused across later capabilities, the tag silently became stale and misleading.

**Resolution:**
Use stable, capability-neutral tags such as `[RAD][RadiationDiagnostics::append]` in shared infrastructure, while reserving historical phase labels for documentation or comments only.

**Prevention Rule:**
Search shared runtime code for hardcoded phase identifiers whenever a diagnostic or helper is reused.

**Lesson:**
Versioned prose may age gracefully; versioned runtime strings in shared code do not.

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
- `bx.smallEnd(2)` and `bx.bigEnd(2)` are the true bounds for the current box
- AMR levels and domains can differ
- Hardcoded bounds create silent correctness bugs

**Common Mistake:**
```cpp
for (int k = 0; k < 50; ++k) {
    tau_cum += tau_per_layer;
}
```

**Lesson:**
Always derive loop bounds dynamically.

---

### B.2 – Vertical Sweep: One Thread per `(i,j)`, Sequential `k` Loop

**Pattern:**
```cpp
const auto& lo = bx.loVect();
const auto& hi = bx.hiVect();
amrex::ParallelFor(amrex::Box(lo[0], lo[1], 0, hi[0], hi[1], 0),
    [=] AMREX_GPU_DEVICE (int i, int j, int) {
        amrex::Real tau_cum = 0.0;
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            tau_cum += tau_per_layer;
        }
    }
);
```

**Why:**
- Horizontal columns are naturally parallel
- Vertical integration is stateful and therefore sequential within a column

**Lesson:**
Map the kernel structure to the underlying mathematical structure.

---

### B.3 – Extract Domain Extents Explicitly

**Pattern:**
```cpp
const auto& domain = geom[lev].Domain();
int kmin = domain.smallEnd(2);
int kmax = domain.bigEnd(2);
```

**Lesson:**
Geometry queries keep the radiation sweep robust across resolution changes and AMR configurations.

---

### B.4 – Preserve a Uniform `dz` Fallback While Preparing for Nonuniform Geometry

**Issue:**
Heating-rate divergence needs physically correct `dz`, but the code must remain safe on legacy uniform grids.

**Resolution:**
Initialize per-level `dz` arrays from the uniform cell size, then override them only when physically resolved vertical geometry is available.

**Lesson:**
Future-ready infrastructure should not break the present default path.

---

## Part C: Atmospheric State Access and Defensive Validation

### C.1 – Safe MultiFab Component Access

**Pattern:**
```cpp
constexpr int Rho_comp = 0;
constexpr int RhoTheta_comp = 4;
constexpr int Temp_comp = 5;

amrex::Real rho = state(i, j, k, Rho_comp);
amrex::Real T   = state(i, j, k, Temp_comp);
if (rho <= 0.0) rho = 1.0;
if (T <= 0.0)   T   = 288.15;
```

**Why:**
- Named indices improve readability and reduce mistakes
- Defensive checks prevent NaN/Inf propagation from bad intermediate states

**Lesson:**
State access should always be explicit, validated, and readable.

---

### C.2 – Temperature and Thermodynamic State Need Explicit Interpretation

**Issue:**
ERF often stores `RhoTheta` rather than temperature directly.

**Resolution:**
When diagnostics or physics require temperature, derive it deliberately from the available state, using the appropriate density and pressure assumptions for the context.

**Lesson:**
Do not assume standard atmospheric variables are stored in their most convenient form.

---

### C.3 – Defensive Clipping for Unphysical Values

**Pattern:**
```cpp
if (T <= 0.0 || T > 400.0) {
    T = 288.15;
}
if (rho <= 0.0) {
    rho = 1.0;
}
```

**Why:**
- Simulations can transiently produce invalid inputs
- Controlled clipping is preferable to silent NaN propagation

**Lesson:**
Defensive physics kernels should fail safe, not fail silently.

---

### C.4 – Logging Without Spamming

**Pattern:**
```cpp
if (verbosity >= 1 && amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Radiation: clipped " << n_clipped << " density values
";
}
```

**Lesson:**
Aggregate errors and report summaries, not per-cell details.

---

### C.5 – Finite Guards Protect Downstream Coupling

**Issue:**
A single NaN in `qheating_rates` can silently corrupt downstream thermodynamic or PBL tendencies.

**Resolution:**
Guard inputs such as `dz`, `rho`, and `cp`; guard outputs such as computed heating rates; and return a conservative safe value when invalid input would otherwise propagate.

**Lesson:**
Finite checks are not optional once radiation outputs feed other physics modules.

---

## Part D: Integration Pitfalls and Resolution Patterns

### D.1 – Pitfall: Hardcoded Grid Bounds

**Issue:**
Early prototypes used fixed vertical limits such as `k < 50`, which becomes unsafe whenever the domain uses a different number of levels.

**Resolution:**
All vertical loops must derive limits from `Box` or geometry metadata.

**Lesson:**
Hardcoded grid bounds are a silent correctness failure.

---

### D.2 – Pitfall: Host Loop Calling a Device Function

**Issue:**
A host-side nested loop once called the device-oriented vertical sweep directly and then accumulated results on host, which is invalid on GPU builds and race-prone on threaded CPU paths.

**Resolution:**
Launch the sweep through `amrex::ParallelFor`, keep reductions on device, and copy back only reduced scalars.

**Lesson:**
Device functions must be invoked only from device execution contexts.

---

### D.3 – Pitfall: Physically Correct Code That Is Not Wired Into the Model

**Issue:**
A radiation driver can be mathematically correct, tested in isolation, and still be inert if it is never called, if its output MultiFab is never allocated, or if downstream consumers ignore its output.

**Resolution:**
Audit four gates together: allocation, driver invocation, per-level output population, and downstream source consumption.

**Lesson:**
Correctness and wiring are separate completion criteria.

---

### D.4 – Pitfall: Temporal Wiring Errors Masquerade as Physics Errors

**Issue:**
If radiation is called too often, too rarely, or at mismatched points relative to source-term application, the symptoms appear as duplicate diagnostics, missing forcing, or unstable tendencies rather than compiler failures.

**Resolution:**
Trace call frequency, source-term gates, and diagnostic cadence together whenever the time-integration path changes.

**Lesson:**
Temporal semantics are part of the physics contract.

---

### D.5 – Pitfall: Diagnostics Identity Weaker Than `(step, time, call_site)`

**Issue:**
Multi-call-site diagnostics cannot safely deduplicate records if identity is keyed on `step` alone.

**Resolution:**
Treat `(step, time, call_site)` as the event identity in the writer, duplicate guard, and validation scripts.

**Lesson:**
Observability infrastructure must evolve with timestep semantics.

---

### D.6 – Pitfall: Inferring Call Site from Out-of-Scope Time Variables

**Issue:**
Lower-level diagnostics code once attempted to infer `call_site` using variables it did not own.

**Resolution:**
Thread semantic labels explicitly from the caller.

**Lesson:**
Semantic context belongs to the call site, not to heuristic reconstruction inside utilities.

---

### D.7 – Pitfall: Validation Scripts Ignore the Configured Diagnostics File

**Issue:**
A checker once searched for a default filename even though the run wrote a case-specific diagnostics file.

**Resolution:**
Validation scripts must follow the configured `diag_file` path or documented expected candidates.

**Lesson:**
An I/O contract mismatch can look like a solver failure.

---

### D.8 – Pitfall: Diagnostics Controls Alter the Physics Path

**Issue:**
Diagnostics cadence or call-site filters must never decide whether heating is computed.

**Resolution:**
Compute physics unconditionally; apply diagnostics controls only when writing or filtering output.

**Lesson:**
Observability settings must remain orthogonal to model evolution.

---

### D.9 – Pitfall: Benchmark Validation Hardcodes One Diagnostics Mode

**Issue:**
Benchmark suites fail spuriously if they assume one row count or one call-site pattern regardless of runtime diagnostics mode.

**Resolution:**
Benchmark case definitions should carry their expected diagnostics cadence, and validation logic must derive row counts and call-site checks from that case metadata.

**Lesson:**
Reproducibility requires configuration-aware validation, not brittle hardcoding.

---

### D.10 – Pitfall: Surface and SEB Ownership Rules Are Implicit

**Issue:**
Surface-property fallback work and later SEB extensions showed that ownership of `t_sfc`, `q_sfc`, and related fields becomes ambiguous when both fallback filling and prognostic updates are active.

**Resolution:**
Document a strict precedence rule: LSM-provided fields when available, otherwise scalar fallbacks, while prognostic mode owns `t_sfc` and `q_sfc` once enabled.

**Lesson:**
Surface-field ownership must be treated as an explicit interface contract.

---

### D.11 – Pitfall: Shared Capability Checks Ignore the TwoStream Path

**Issue:**
Plotfile capability checks once recognized only the non-TwoStream radiation path when deciding whether radiation heating output existed.

**Resolution:**
Capability checks must use the same gating logic as `qheating_rates` allocation so that advertised outputs match actual storage.

**Lesson:**
Shared introspection utilities must track all active physics branches.

---

### D.12 – Pitfall: Cloud-Fraction Diagnosis Without Physical Bounds

**Issue:**
RH- and `qc`-based cloud fraction can become numerically unstable or unphysical if evaluated without bounds and finite guards.

**Resolution:**
Blend RH and condensate indicators, clamp cloud fraction to `[0,1]`, and protect all saturation calculations against invalid thermodynamic states.

**Lesson:**
Cheap diagnostic closures still require physical and numerical guardrails.

---

## Part E: Testing and Validation Checklist

### E.1 – Compilation and GPU Readiness
- [ ] Code compiles with the intended GPU toolchain and architecture settings
- [ ] No warnings indicate invalid host-device execution patterns
- [ ] No warnings indicate unintended host-device copies or uninitialized device memory

### E.2 – Grid Adaptivity
- [ ] Coarse and fine grids both produce valid diagnostics
- [ ] Vertical loop counts match the actual grid layout
- [ ] No hardcoded bounds remain in the active code path

### E.3 – Physical Correctness
- [ ] SW flux decreases monotonically with increasing optical depth where expected
- [ ] LW fluxes remain smooth and finite for valid thermodynamic states
- [ ] Heating rates have physically reasonable sign and magnitude
- [ ] No NaN/Inf appears in diagnostics or coupled source fields
- [ ] Scattering-only controls behave consistently when scattering albedo is zero versus nonzero

### E.4 – Wiring and Coupling
- [ ] Driver call sites are confirmed by code search, not inferred from local edits alone
- [ ] Output MultiFabs are allocated under the same conditions that gate the driver call
- [ ] Downstream consumers read the same outputs under matching gates

### E.5 – Diagnostics Cadence and Identity
- [ ] Diagnostics schema includes `call_site` when multiple call sites are supported
- [ ] Duplicate suppression keys on `(step, time, call_site)`
- [ ] Validation scripts respect configured diagnostics mode and diagnostics file path
- [ ] Case expectations for rows per step are explicit and tested

### E.6 – Regression and Reproducibility
- [ ] Shared runtime code is checked for stale phase identifiers
- [ ] Existing canonical cases remain unchanged unless the new capability is explicitly enabled
- [ ] Benchmark tolerances are centralized and version controlled
- [ ] Machine-readable and human-readable reports remain deterministic

---

## Part F: Documentation Standards and Summary Checklist

### F.1 – Function-Level Documentation

Document GPU-safe helper functions with equations, units, bounds, and fallback behavior so that reviewers can verify both physics intent and device suitability.

### F.2 – File Headers

File headers should describe the supported capability set, the governing approximations, and the primary references without embedding transient task-tracking language.

### F.3 – Code Comments

Code comments should explain numerical intent, ordering constraints, or interface contracts. Historical notes are acceptable in comments and design documents, but runtime behavior should remain phase-agnostic.

### F.4 – Key Principles

1. Mark device-callable helpers explicitly.
2. Keep host I/O out of device kernels.
3. Derive loop bounds from runtime geometry.
4. Match kernel structure to column-physics structure.
5. Validate and clip unsafe atmospheric inputs conservatively.
6. Keep shared runtime strings capability-neutral.
7. Verify wiring, allocation, and downstream consumption together.
8. Treat diagnostic identity as `(step, time, call_site)` when cadence is multi-site.
9. Keep diagnostics controls observational, never dynamical.
10. Make validation scripts configuration-aware and reproducible.

---

## References

- AMReX GPU Guide: https://amrex-codes.github.io/amrex/docs_html/GPU.html
- AMReX Parallel Loop Patterns: https://amrex-codes.github.io/amrex/docs_html/GPU_HowTo.html#parallel-for
- Atomic Operations: https://amrex-codes.github.io/amrex/docs_html/GPU_HowTo.html#atomic-operations
- Meador, W. E., and W. R. Weaver, 1980: "Two-stream approximations to radiative transfer in planetary atmospheres", *J. Atmos. Sci.*, 37, 630–643.

---

## Part D Supplement: Diagnostics Dedup and Nonuniform `dz` Lessons

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

### D.13 – Summary of Diagnostics Dedup and Nonuniform `dz` Guidance

**Dedup + Nonuniform dz + Finite checks** = robust integration polish that:
- Prevents silent data corruption (dedup tuples, finite checks)
- Prepares for terrain-aware grids without breaking uniform-grid behavior
- Remains GPU-safe and performance-neutral on current (uniform) simulations

---

## Part D: Phase 14 Lessons (Prognostic Cloud Fraction)

### D.1 – RH/qc-based Cloud Fraction Diagnosis: Physical Consistency

**Pattern (Phase 14):**
```cpp
// Diagnose RH from water vapor mixing ratio and T
amrex::Real rh = compute_relative_humidity(qv, T, P);

// Cloud fraction from RH ramp + qc scaling
amrex::Real cf = diagnose_cloud_fraction_from_rh_qc(
    rh, qc, rh_min, rh_max, qc_scale);

// Use cf(k) to scale cloud optical depth at each level
tau_cloud(k) = cf(k) * cloud_tau_per_layer;
```

**Why:**
- RH provides a thermodynamic signal (saturation proximity); qc provides microphysical evidence
- Per-level cf(k) is more physical than global scalar cloud_fraction
- Scaling (not binary on/off) allows smooth transitions and sub-grid variability
- Linear RH ramp + qc blend is computationally efficient and tunable

**Common Mistake:**
```cpp
// ❌ WRONG: Hardcoded cloud fraction threshold
if (qc > 1e-5) cf = 1.0;  // Binary; no gradation with RH

// ❌ WRONG: Unbounded cf from qc alone
cf = qc_scale * qc;  // Can exceed 1.0 if qc_scale too large

// ❌ WRONG: RH without bounds or fallback
rh = qv / qsat;  // Divide by zero if qsat=0 or invalid T
```

**Lesson:**
- Blend RH and qc signals for physical robustness
- Always clamp cf ∈ [0, 1] and guard against NaN/Inf
- Test sensitivity to rh_min, rh_max, qc_scale; provide safe defaults

---

### D.2 – Finite Guards in RH Computation: Magnus Saturation Formula Safety

**Pattern (Phase 14):**
```cpp
// Saturation vapor pressure (Magnus formula)
amrex::Real e_sat = e0 * std::exp(a * (T - T0) / (T - b));
amrex::Real qsat = epsilon * e_sat / (P - e_sat);

// Guard against overflow and division by zero
if (arg > 100.0) arg = 100.0;  // Prevent std::exp overflow
if (qsat < 0.0 || !std::isfinite(qsat)) qsat = 1.0e-6;  // Fallback

// RH clamped to [0, 1]
amrex::Real rh = amrex::min(1.0, amrex::max(0.0, qv / qsat));
```

**Why:**
- Magnus exp argument can overflow if T >> reference (defensive clipping needed)
- qsat → ∞ as P → e_sat (division by zero risk near saturation)
- qv/qsat can exceed 1.0 or become NaN if qsat invalid (must clamp)
- RH approximations work only in reasonable temperature range (~200–350 K)

**Common Mistake:**
```cpp
// ❌ WRONG: No clipping on exp argument
amrex::Real e_sat = e0 * std::exp(a * (T - T0) / (T - b));  // Overflow if T very large

// ❌ WRONG: No fallback for qsat ≈ 0
amrex::Real rh = qv / qsat;  // NaN if qsat=0 or very small
```

**Lesson:**
- Clip exp arguments to avoid overflow; fallback to a large value (e.g., 100)
- Test with T outside normal range (e.g., T = 500 K or T = 100 K) to verify guards
- Ensure qsat > 0 before division; use reasonable fallback (e.g., 1 kg/kg)

---

### D.3 – Per-Level Diagnosis in Vertical Sweeps: Integration Complexity

**Pattern (Phase 14):**
```cpp
// Inside vertical sweep loop (SW down, LW up/down):
for (int k = kmin; k <= kmax; ++k) {
    // Compute state-dependent optical depth for THIS level
    amrex::Real tau = tau_layer_value(k, ...);  // Base tau
    
    if (rad_choice.tau_sw_dynamic_enable) {
        tau = diagnose_tau_sw_dynamic(i, j, k, state_arr, tau, rad_choice);
    }
    
    // NEW Phase 14: Per-level cloud fraction modulation
    if (rad_choice.cloud_fraction_prog_enable && cloudy) {
        amrex::Real cf = diagnose_cloud_fraction_prognostic(i, j, k, state_arr, rad_choice, geom);
        amrex::Real tau_base = tau_layer_value(k, ..., /*cloudy=*/false);
        tau = tau_base + cf * cloud_tau_per_layer;  // Scale cloud component
    }
    
    // Apply tau to flux computation
    amrex::Real flux = compute_flux(tau, ...);
    // ... accumulate ...
}
```

**Why:**
- Each level has different qv, qc, T → different RH and thus cf
- Scaling cloud_tau_per_layer by cf(k) makes cloud impact physically local
- Integration point must be after other tau diagnostics (Phase 12 dynamic tau)
- Order matters: Phase 12 first, then Phase 14, then flux computation

**Common Mistake:**
```cpp
// ❌ WRONG: Compute cf once for whole column, apply globally
amrex::Real cf_global = diagnose_cf(...);  // Only once, outside loop
for (int k = kmin; k <= kmax; ++k) {
    tau = tau_base + cf_global * cloud_tau_per_layer;  // Wrong: ignores level differences
}

// ❌ WRONG: Apply cf before dynamic tau, losing dynamic effect
tau = tau_base + cf * cloud_tau_per_layer;  // Phase 14
tau = diagnose_tau_sw_dynamic(..., tau, ...);  // Phase 12 overwrites!
```

**Lesson:**
- Diagnose cf(k) at EVERY level; never hoist out of sweep loop
- Order integration steps: Phase 3 (base) → Phase 12 (dynamic) → Phase 14 (prognostic cf)
- Test with varying qv/qc profiles (e.g., dry below 500m, cloud layer 500–2000m) to verify per-level diagnosis

---

### D.4 – Temporal Smoothing State: Future Infrastructure Requirement

**Issue (Phase 14):**
```cpp
// Current implementation (Phase 14): NO temporal smoothing at sweep level
if (rad_choice.cloud_fraction_smooth_enable && rad_choice.cloud_fraction_smooth_alpha > 0.0) {
    // ❌ BLOCKED: No persistent storage for cf_old(i, j, k) across timesteps
    // Would require new MultiFab (like qheating_rates) initialized in ERF_MakeNewArrays.cpp
}

// Smoothing logic is ready (in ERF_PrognosticCloudFraction.H):
cf_smooth = smooth_cloud_fraction_ema(cf_new, cf_old, alpha);
// But cf_old is not available in vertical_two_stream_sweep()
```

**Why:**
- EMA smoothing requires storing cf from previous timestep: cf_smooth(i,j,k,t-dt)
- vertical_two_stream_sweep() is device-side with no access to persistent storage
- Adding per-level state MultiFab requires coordination with ERF initialization and boundary handling
- Phase 14 focuses on diagnosis; smoothing deferred to Phase 15+ infrastructure work

**Prevention:**
- Set `cloud_fraction_smooth_enable = false` (default) to disable smoothing
- Set `cloud_fraction_smooth_alpha = 0.0` (default) to disable smoothing
- Parameters validated but NOT applied; future phases can add state storage

**Lesson:**
- Persistent per-level state in radiation solvers requires MultiFab infrastructure
- Device-side code cannot dynamically allocate or access time-history arrays
- Plan MultiFab addition (ERF_MakeNewArrays.cpp, boundary handling, diagnostics export) in separate phase
- Current Phase 14 validates parameters and keeps smoothing logic ready for future use

---

### D.5 – Finite Guards: Saturation and Defensive Fallbacks

**Pattern (Phase 14):**
```cpp
// Diagnose cloud fraction with multiple defensive layers
amrex::Real cf = diagnose_cloud_fraction_from_rh_qc(rh, qc, rh_min, rh_max, qc_scale);

// Inside diagnose_cloud_fraction_from_rh_qc():
// Guard 1: Input validation
if (!std::isfinite(rh) || !std::isfinite(qc)) return 0.0;

// Guard 2: Physical range checks
if (rh < 0.0) rh = 0.0;
if (rh > 1.0) rh = 1.0;
if (qc < 0.0) qc = 0.0;

// Guard 3: Compute with saturation (min(1.0, x))
amrex::Real cf_rh = (rh - rh_min) / (rh_max - rh_min);  // Linear ramp
if (cf_rh < 0.0) cf_rh = 0.0;
if (cf_rh > 1.0) cf_rh = 1.0;

amrex::Real cf_qc = qc_scale * qc;
if (cf_qc > 1.0) cf_qc = 1.0;

amrex::Real cf = cf_rh + cf_qc;  // Sum saturates
if (cf > 1.0) cf = 1.0;  // Guard 4: Final clamp

// Guard 5: Finiteness check before return
if (!std::isfinite(cf)) cf = 0.0;
return cf;
```

**Why:**
- qv/qc/T from state arrays may be uninitialized, NaN, or corrupted in edge cases
- RH formula can fail if P ≈ e_sat (saturation); clamping prevents division singularities
- cf ∈ [0, 1] is hardened at multiple levels (input, intermediate, output)
- Silent fallback to 0 (no cloud fraction) is conservative; simulation continues

**Common Mistake:**
```cpp
// ❌ WRONG: Minimal guards
amrex::Real cf = (rh - rh_min) / (rh_max - rh_min) + qc_scale * qc;
// Missing: divide by zero if rh_max==rh_min, no bounds on cf, no NaN check

// ❌ WRONG: Late clamping only
amrex::Real cf = rh_contribution + qc_contribution;
if (cf > 1.0) cf = 1.0;  // Too late; cf > 1 may have already corrupted downstream code
```

**Lesson:**
- Apply guards at input, intermediate steps, and output (defense in depth)
- Always check `std::isfinite()` for any quantity derived from state arrays
- Saturate (clamp) immediately after combining contributions (RH + qc)
- Test with invalid states (T=0, P<0, qv=NaN, qc=Inf) to verify guards catch all paths

---

### D.15 – Summary of Prognostic Cloud-Fraction Guidance

**Prognostic Cloud Fraction** combines thermodynamic (RH) and microphysical (qc) signals with:
- Finite-guarded RH diagnosis from Magnus saturation formula
- Per-level cf(k) diagnosis in vertical sweeps (after Phase 12 dynamic tau)
- Saturation blending (cf ≤ 1) of RH + qc contributions
- Temporal smoothing infrastructure ready for Phase 15+ persistent-state addition
- Multiple defensive guards at input, intermediate, and output stages

Enables physically consistent, per-level cloud fraction modulation of radiation while maintaining backward compatibility and GPU safety.

---

## Part D.3 – Phase 14B Lesson: Multilevel Vector Allocation for Radiation State

### D.3.1 – Allocate Before Use: Radiation Surface-Property MultiFabs

**Pattern:**
```cpp
// ERF_Constructors.cpp::ERF_shared() – Constructor stage
int nlevs_max = max_level + 1;

// Step 1: Resize all radiation state vectors
qheating_rates.resize(nlevs_max);
rad_fluxes.resize(nlevs_max);
twostream_alb_sw.resize(nlevs_max);      // Phase 14A/14B: Fallback surface albedo
twostream_emiss_lw.resize(nlevs_max);    // Phase 14A/14B: Fallback surface emissivity
twostream_t_sfc.resize(nlevs_max);       // Phase 14A/14B: Fallback surface temperature

// Later in init_stuff() – Per-level allocation stage
if (solverChoice.radChoice.rad_type == RadType::TwoStream) {
    twostream_alb_sw[lev]   = std::make_unique<MultiFab>(ba2d[lev], dm, 1, 0);
    twostream_emiss_lw[lev] = std::make_unique<MultiFab>(ba2d[lev], dm, 1, 0);
    twostream_t_sfc[lev]    = std::make_unique<MultiFab>(ba2d[lev], dm, 1, 0);
    
    // Initialize from scalar fallbacks
    twostream_alb_sw[lev]->setVal(solverChoice.radChoice.surface_albedo_sw);
    twostream_emiss_lw[lev]->setVal(solverChoice.radChoice.surface_emissivity_lw);
    twostream_t_sfc[lev]->setVal(solverChoice.radChoice.surface_temp_k);
}
```

**Why:**
- `Vector<std::unique_ptr<MultiFab>>` must be **pre-sized in constructor** (ERF_shared)
- **Before** entering per-level allocation loop (init_stuff)
- If `.resize(nlevs_max)` is skipped, indexing `vector[lev]` in init_stuff is **undefined behavior**:
  - Vector capacity is 0; accessing index [lev] accesses unallocated memory
  - May crash, silently corrupt state, or hang on GPU
  - Bug not caught at compile-time (Vector allows operator[] without bounds checking)

**Common Mistake:**
```cpp
// ❌ WRONG: Forgot to resize in constructor
// ERF_Constructors.cpp::ERF_shared()
// (no resize call for twostream_alb_sw)

// Later in init_stuff():
if (solverChoice.radChoice.rad_type == RadType::TwoStream) {
    // CRASH: twostream_alb_sw.size() == 0, accessing [lev] is undefined behavior!
    twostream_alb_sw[lev] = std::make_unique<MultiFab>(ba2d[lev], dm, 1, 0);
}
```

**Fix Checklist:**
1. Identify all radiation state vectors that hold per-level data (MultiFab, FAB, Array, etc.)
2. For each vector, add a `.resize(nlevs_max)` call in `ERF::ERF_shared()` **after** `int nlevs_max = max_level + 1`
3. Place all radiation vector resizes together (near qheating_rates/rad_fluxes) for maintainability
4. Comment with the phase that added the vector (e.g., `// Phase 14A/14B`)
5. Verify in init_stuff() that vector indices are now valid before calling `.make_unique<>()`

### D.3.2 – Conditional Allocation Requires Unconditional Sizing

**Pattern:**
```cpp
// Size unconditionally in constructor (regardless of whether allocation will happen)
twostream_alb_sw.resize(nlevs_max);  // Done even if TwoStream is disabled!

// Later, allocate conditionally per-level
for (int lev = 0; lev <= max_level; ++lev) {
    init_stuff(lev, ...);  // Inside: if (radChoice.rad_type == TwoStream) { allocate }
}
```

**Why:**
- Sizing is cheap (just reserves capacity)
- Allocation (`.make_unique<>()`) can be expensive and conditional on runtime parameters
- Decoupling them prevents hard-to-debug indexing bugs
- Future phase may enable TwoStream conditionally; size() call is already there

**Backward Compatibility:**
- If TwoStream is disabled, vectors are sized but empty (unique_ptrs are nullptr)
- No performance cost (empty vectors are negligible)
- No numerical impact (code never accesses nullptr vectors when TwoStream is off)

### D.3.3 – Phase 14B Regression Test Pattern

**Test Setup** (TwoStream_ProgCloudFraction):
```
# Test inputs file: erf.radiation.rad_type = TwoStream
# Surface properties: fallback scalars (no LSM)
# Grid: Simple 32×32×16 box; 5 timesteps
# Purpose: Verify TwoStream vectors are properly sized and allocated
```

**Validation Checks:**
1. **No crash on startup**: Constructor properly sizes vectors
2. **Output bitwise-identical to baseline**: Surface properties via fallback MultiFabs match pre-14A scalar-only runs
3. **Quiet operation**: No NaN warnings or invalid-value messages

---

### D.16 – Summary of Surface-Property Wiring Guidance

**Multilevel Vector Allocation** for radiation state must follow the two-stage pattern:
- **Stage 1 (Constructor)**: Resize all vectors in `ERF::ERF_shared()` unconditionally
- **Stage 2 (Init)**: Allocate individual level MultiFabs conditionally in `init_stuff()` per-level
- **Decoupling** avoids indexing bugs while maintaining flexibility for conditional features
- **Comment and document** why each radiation vector needs resizing (for future maintainers)

Enables safe, scalable per-level radiation state management across multiple AMR levels.

---

## Part D Supplement: Aerosol Optical-Depth Integration Lessons

**Status**: ✅ Complete (as of 2026-08-09)  
**Key Skills Demonstrated**: Device-safe profile-based parameter handling, enum-based runtime selection, height-dependent kernel logic

### Lesson: Height-Dependent Kernel Logic Without External Array Access

**Pattern:**
```cpp
// In vertical_two_stream_sweep():
for (int k = kmin; k <= kmax; ++k) {
    // Compute height incrementally
    amrex::Real z_level = 0.0;
    for (int kk = kmin; kk < k; ++kk) {
        z_level += dz_uniform;
    }
    
    // Call device function with height parameter
    if (rad_choice.aerosol_profile_type == AerosolProfileType::Exponential) {
        tau_aerosol = diagnose_tau_aerosol_exponential(z_level, tau_surface, scale_height);
    }
}

// Device-side (ERF_AerosolOpticalDepth.H):
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_tau_aerosol_exponential(amrex::Real z_level,
                                            amrex::Real tau_surface,
                                            amrex::Real scale_height_m)
{
    // Guard against invalid inputs
    if (!std::isfinite(z_level) || !std::isfinite(tau_surface) || !std::isfinite(scale_height_m)) {
        return 0.0;
    }
    if (scale_height_m <= 0.0 || tau_surface < 0.0) {
        return 0.0;
    }
    
    // Compute with overflow protection
    amrex::Real arg = -z_level / scale_height_m;
    if (arg < -100.0) return 0.0;  // exp would underflow
    
    amrex::Real tau = tau_surface * std::exp(arg);
    
    // Clamp output to [0, 100]
    if (tau < 0.0) tau = 0.0;
    if (tau > 100.0) tau = 100.0;
    
    return tau;
}
```

**Why This Approach:**
1. **No external array access** inside device function; all parameters passed by value
2. **Overflow protection**: exp() argument guarded to avoid inf/nan
3. **Finite guards**: All inputs validated before computation
4. **Clamped output**: Result guaranteed in [0, 100] physical range
5. **GPU-safe**: Pure computation, no host I/O or dynamic allocation

**Common Mistake:**
```cpp
// ❌ WRONG: Pass array, access inside device function
for (int k = kmin; k <= kmax; ++k) {
    tau_aerosol = diagnose_tau_aerosol(state_arr, k);  // Might fail on GPU!
}

// Later, in device function:
AMREX_GPU_DEVICE
amrex::Real diagnose_tau_aerosol(const Array4<amrex::Real>& state_arr, int k) {
    // ❌ WRONG: Array bounds not checked in device context
    amrex::Real z = state_arr(i, j, k, zcomp);  // i, j undefined!
}
```

**Lesson:** Precompute derived values (like height) in host loop, pass as scalar parameters to device function. Avoid passing arrays and accessing them inside device kernels when simpler scalar approach suffices.

### Lesson: Enum-Based Runtime Profile Selection on GPU

**Pattern:**
```cpp
// Host-side: ERF_RadStruct.H
AMREX_ENUM(AerosolProfileType, Constant, Exponential, Table);

// In RadChoice struct:
AerosolProfileType aerosol_profile_type = AerosolProfileType::Constant;

// Device-side: vertical_two_stream_sweep()
for (int k = kmin; k <= kmax; ++k) {
    amrex::Real tau_aerosol = 0.0;
    
    if (rad_choice.aerosol_profile_type == AerosolProfileType::Constant) {
        tau_aerosol = diagnose_tau_aerosol_constant(rad_choice.aerosol_tau_per_layer);
    } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Exponential) {
        amrex::Real z_level = ...;  // Compute height
        tau_aerosol = diagnose_tau_aerosol_exponential(z_level, rad_choice.aerosol_tau_surface, 
                                                      rad_choice.aerosol_scale_height_m);
    } else if (rad_choice.aerosol_profile_type == AerosolProfileType::Table) {
        tau_aerosol = diagnose_tau_aerosol_table(k);
    }
    
    tau_sw += tau_aerosol;
}
```

**Why This Pattern:**
1. **AMREX_ENUM macro**: Generates both host and device versions safely
2. **Branch in host loop**: Profile type checked once per level in outer (host-side) loop
3. **Device functions branchless**: Each device function focuses on single profile, no conditional
4. **Future extensibility**: Adding new profile type requires only new function + new enum value + new branch
5. **GPU performance**: Host-loop branching (cheap) vs device branching (expensive in divergent warp) avoided

**Common Mistake:**
```cpp
// ❌ WRONG: String comparison in device kernel
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_tau_aerosol_gpu(const char* profile_type, ...) {
    if (std::string(profile_type) == "exponential") {  // ❌ std::string not device-safe!
        ...
    }
}
```

**Lesson:** Use AMREX_ENUM for device-safe type selection. Branch on enum in host loop, not in device kernel.

### Lesson: Additive Physics Integration Without Double-Counting

**Pattern:**
```cpp
// Compute base optical depth (from Phase 1-11)
amrex::Real tau = tau_layer_value(...);

// Add Phase 12 dynamic tau (if enabled)
if (rad_choice.tau_sw_dynamic_enable) {
    tau = diagnose_tau_sw_dynamic(...);  // Returns tau_base + dynamic_component
}

// Add Phase 14 cloud fraction scaling (if enabled and cloudy)
if (rad_choice.cloud_fraction_prog_enable && cloudy && ...) {
    amrex::Real cf = diagnose_cloud_fraction_prognostic(...);
    amrex::Real tau_base_k = tau_layer_value(..., /*cloudy=*/false);
    tau = tau_base_k + cf * rad_choice.cloud_tau_per_layer;
}

// Add Phase 15 aerosol (if enabled)
if (rad_choice.aerosol_enable) {
    amrex::Real tau_aero = diagnose_tau_aerosol_*(...);
    tau += tau_aero;  // ✅ ADDITIVE: += not =
}

// Use combined tau in flux computation
F_sw = compute_sw_flux(tau, ...);
```

**Why Additive:**
1. **Separation of concerns**: Each phase contributes independently
2. **No double-counting**: Cloud fraction already accounted for in tau_base_k recovery above
3. **Backward compat**: When Phase 15 disabled, aerosol contribution is exactly zero
4. **Physical interpretation**: tau_total = clear-sky + cloud + dynamic moisture + aerosol

**Common Mistake:**
```cpp
// ❌ WRONG: Replacing instead of adding
if (rad_choice.aerosol_enable) {
    tau = diagnose_tau_aerosol(...);  // ❌ Replaces, loses cloud/dynamic contributions!
}
```

**Lesson:** When integrating new physics contributions, use `+=` to stack on top, not `=` to replace.

### Lesson: Parameter Validation and Default Safeguards

**Pattern:**
```cpp
// In init_params() (ERF_RadStruct.H):
pp.query("radiation.aerosol_tau_per_layer", aerosol_tau_per_layer);
pp.query("radiation.aerosol_scale_height_m", aerosol_scale_height_m);
pp.query("radiation.aerosol_tau_surface", aerosol_tau_surface);

// Validation: clamp to physically reasonable ranges
if (aerosol_tau_per_layer < 0.0) aerosol_tau_per_layer = 0.0;
if (aerosol_tau_surface < 0.0) aerosol_tau_surface = 0.0;
if (aerosol_scale_height_m <= 0.0) aerosol_scale_height_m = 2000.0;  // Restore default
```

**Why This Approach:**
1. **Soft clipping**: Invalid inputs don't crash; they're silently corrected to safe defaults
2. **Defensive coding**: User mistakes (negative tau, zero scale height) don't propagate into kernels
3. **Device-side simplification**: Device functions assume inputs already validated
4. **Diagnostic output**: On host after init, invalid params are logged/corrected before simulation

**Testing Pattern:**
```python
# In RegTest check script:
for fval in [tau, flux, heating_rate]:
    if fval != fval:  # NaN check
        print(f"ERROR: NaN in diagnostics")
        return False
    if abs(fval) == float('inf'):  # Inf check
        print(f"ERROR: Inf in diagnostics")
        return False
```

**Lesson:** Validate all user input in init_params(). Device functions should assume input is already clean.

### D.14 – Summary of Aerosol Optical-Depth Guidance

**Profile-Based Physics Stacking** for complex multi-component models:
1. Define profile type via AMREX_ENUM (host+device safe)
2. Precompute scalar parameters in host loop (e.g., height)
3. Call profile-specific device functions with scalar params (no array access)
4. Stack contributions additively with `+=` (not replacement)
5. Clamp/validate all parameters in init_params(); device functions assume clean input
6. Guard all device logic against NaN/Inf/invalid bounds with early returns

Enables modular, extensible, GPU-safe integration of multi-phase physics without sacrificing clarity or safety.
