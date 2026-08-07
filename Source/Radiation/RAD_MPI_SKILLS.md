Here's the `RAD_MPI_SKILLS.md` content again, cleanly, for manual copy-paste:

```markdown name=RAD_MPI_SKILLS.md
# ERF Two-Stream Radiation Module — Development Skills & Bug Fix Reference

Complete record of lessons learned during development of the ERF Two-Stream
Radiation module. Use this as a checklist before merging any new Radiation
phase. Modeled directly on `Source/UrbanCanopy/UCM_MPI_SKILLS.md` — same
categories, same discipline, applied to the Radiation module's specific
bugs and conventions.

This file will grow phase by phase as new technical challenges are
discovered and documented.

---

## Part A — Architecture & Design Rules

### A1. Follow the UCM/Dust Module Wiring Patterns Exactly

Every new Radiation source file must be registered in **both** build
systems, exactly as required for UCM/Dust (see `UCM_MPI_SKILLS.md` A2).

**Make.package pattern:**

```
CEXE_sources += Radiation/ERF_NewRadiationFile.cpp
CEXE_headers += Radiation/ERF_NewRadiationFile.H
```

**CMake pattern (in `CMake/BuildERFExe.cmake`):**

```
target_sources(${erf_lib_name} PRIVATE ${SRC_DIR}/Radiation/ERF_NewRadiationFile.cpp)
```

**Lesson (Phase 1b):** A fully correct file
(`Source/DataStructs/ERF_RadStruct.H`) sat completely unused for an entire
merged PR because nothing included it and nothing called
`init_params(...)`. Neither CMake registration nor `SolverChoice` wiring
were checked before merge. **Rule:** every PR must grep-confirm both:

```
grep -rn "radChoice" Source/DataStructs/ERF_DataStruct.H
grep -rn "ERF_RadiationDiagnostics.cpp" CMake/BuildERFExe.cmake Source/Radiation/Make.package
```

and paste the output in the PR description (see Contract R2 in
`RAD_DEVELOPMENT.md`).

### A2. Radiation Type Is Opt-In — Never Break the None Default

`RadType::None` must remain the default, and all `erf.radiation.*`
sub-option ParmParse queries must be gated behind
`if (rad_type == RadType::TwoStream)`. This mirrors UCM Contract #37
("optionality: opt-in preservation") — every new Radiation phase must be
byte-identical to the previous phase's output when radiation is disabled,
and any existing non-Radiation ERF test must be completely unaffected.

**Grep check:**

```
grep -n "pp.query(\"radiation." Source/DataStructs/ERF_RadStruct.H
```

Every line found here must be inside the `if (rad_type ==
RadType::TwoStream)` block, never unconditional.

### A3. Diagnostic Stub Call Sites Must Scale With Real Grid State

**Lesson (Phase 1c):** The single most expensive bug in Phase 1 was a
one-line stub (`tau_cum = rad_choice.tau_per_layer;`) that ignored the
actual number of vertical layers entirely. The diagnostic signature of
this bug class is: the wrong answer does not change when you change
the grid resolution. This is the Radiation-module analog of UCM Lesson
18 ("numerical impossibility is a coefficient bug, not a physics bug") —
if a computed flux/quantity is suspicious, first check whether it responds
correctly to resolution changes before hypothesizing about missing
physics.

**Rule:** Any per-column loop or "quantity accumulated over N layers" must
derive N from `geom[lev].Domain().length(2)` (or the equivalent real
box/column extent) — never a hardcoded constant, and never a value that
happens to work for one specific test's grid size.

**Diagnostic pattern to use when debugging a suspicious flux value:**

```
1. Compute the flux/quantity at two different grid resolutions (e.g., n_cell_z = 20 vs 64)
2. If the result is IDENTICAL despite the resolution change, the bug is a
   hardcoded/stub value, not a subtle numerical error
3. If the result changes but not proportionally to the analytical
   expectation, check loop bounds and accumulation logic
```

### A4. RegTest Physical Realism Must Be Checked By Hand, Not Just Numerically

**Lesson (Phase 1, manual fix):** A RegTest passed its own internal
Beer-Lambert self-consistency check while representing a physically
absurd atmosphere (surface SW flux of ~1 W/m² for a "clear-sky" test).
Numerical self-consistency (the solver matches its own analytical formula)
is necessary but not sufficient — the chosen input parameters must
also be sanity-checked against known physical ranges:

| Quantity | Physically plausible range |
|---|---|
| Clear-sky broadband SW optical depth (tau_total) | ~0.1 – 0.3 |
| Clear-sky surface SW flux (moderate zenith angle) | several hundred W/m² |
| Tropospheric LW flux (up or down) | a few hundred W/m² |
| Floating-point relative tolerance (check scripts) | ~1e-6 to 1e-4, never tighter |

**Rule:** After choosing tau_per_layer, tau_lw_per_layer, or any other
per-layer optical property alongside a chosen amr.n_cell z-resolution,
compute tau_total = tau_per_layer * n_layers by hand and confirm the
resulting analytical flux is physically plausible before finalizing the
inputs file. If a domain's vertical resolution is later changed, the
per-layer optical property MUST be re-derived to preserve a sane
tau_total — do not silently let a stale per-layer constant produce
a wildly different total.

### A5. Sounding and Input File Conventions Must Mirror Exec/CanonicalTests/ABL/*

**Lesson (Phase 1, manual fix, twice):**

1. Sounding files must live inside their own case folder, named after
   that specific case (e.g., input_sounding_sw_clearsky inside
   SW_ClearSky_Analytical/) — never shared at the parent module directory
   and referenced by multiple cases via relative path.
2. Sounding file format is fixed: line 1 is a 3-column header
   (Ps/z_ref Ts Qv_s-style surface reference); subsequent lines are
   5 columns in order z [m], theta [K], qv [kg/kg], u [m/s], v [m/s].
   Always fetch a real reference sounding from the repo
   (e.g. Exec/CanonicalTests/ABL/MRF_YSUNew_Enhancements/canonical/sounding_neutral_abl)
   before writing a new one — never assume the column layout.
3. inputs files must match the full structure of an established
   canonical reference (e.g.
   Exec/CanonicalTests/ABL/MRF_YSUNew_Enhancements/canonical/neutral_abl):
   descriptive comment header, physically meaningful stop_time, full
   surface_layer MOST block, PBL parameter block (with commented-out
   alternates as documentation), Coriolis + geostrophic forcing block,
   erf.data_log/erf.profile_int diagnostics — not a stripped-down,
   minimal input file.
4. When moisture is required, enable it the minimal way
   Exec/CanonicalTests/Bomex/inputs_bomex does:
   erf.moisture_model = "SAM" (or a simpler variant if appropriate) +
   erf.buoyancy_type = 1, with a real non-zero qv column in the
   sounding — do not restructure the whole input file.

**Grep check before finalizing any new Radiation RegTest:**

```
# Sounding file must be inside the case folder, not shared at parent level
ls Exec/CanonicalTests/Radiation/<CaseName>/input_sounding_*
# Should NOT find a stray shared sounding file at the parent Radiation/ level
ls Exec/CanonicalTests/Radiation/*.{txt,sounding} 2>/dev/null && echo "WARNING: shared sounding file found outside case folders"
```

---

## Part B — GPU Safety Rules (applies from Phase 2 onward)

### B1. Per-Column Vertical Sweeps Must Be a Single Kernel, Not One Launch Per Level

The two-stream vertical sweep is inherently sequential in k (each level
depends on the cumulative optical depth from levels above/below it).
**Rule:** implement it as one amrex::ParallelFor with one GPU thread per
(i,j) column, looping over the full physical k-range internally within
the device lambda — never launch a separate kernel per vertical level.

```
// CORRECT pattern:
ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int /*k unused*/) noexcept
{
    Real tau_cumulative = 0.0;
    for (int k = kmax; k >= kmin; --k) {
        tau_cumulative += tau_per_layer;   // real accumulation, not a stub
        flux_arr(i,j,k) = compute_sw_direct_flux(tau_cumulative, S0, cos_zenith);
    }
});
```

### B2. No Host-Side Debug Output Inside Device Kernels

Reduce diagnostic scalars (surface flux, TOA flux, max heating rate) on
device first (e.g., via a scratch MultiFab + host-side .max()/.min()
reduction, or amrex::ReduceOps), copy back to host, then print/log from
host only — guarded by amrex::ParallelDescriptor::IOProcessor(). This
mirrors UCM Part B (MPI Collectives Must Come Before IOProcessor Guard) —
any MPI-aware reduction must run on all ranks before any single-rank print.

### B3. Loop Bounds From Real Geometry, Never Hardcoded

Every vertical loop bound must come from geom[lev].Domain() or the
actual box/column extent passed into the kernel — see Contract R3 in
RAD_DEVELOPMENT.md and Part A3 above. This is the single most repeated
lesson across Phase 1/1b/1c and must not recur in Phase 2+.

### B4. Mark All Inline Helpers AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE

Following UCM Lesson 28 (UCM_MPI_SKILLS.md): if a function is intended
to be inlined at every call site (no external symbol), it must be defined
in the header, not split between a .H declaration and a .cpp
definition — mixing the two produces linker errors (one translation unit
sees the inline definition and emits no symbol, another sees only the
declaration and emits a call to an undefined external symbol).

**Grep merge-blocker for future Radiation files:**

```
for cpp in Source/Radiation/*.cpp; do
    if grep -q 'AMREX_FORCE_INLINE' "$cpp"; then
        echo "FAIL: FORCE_INLINE in .cpp file: $cpp — move to header"
    fi
done
```

---

## Part C — Testing and Validation

### C1. Always Compute the Analytical Reference By Hand Before Finalizing a RegTest

Before finalizing any inputs file, manually compute:

1. The expected analytical flux/quantity given the chosen parameters.
2. Confirm it is physically plausible (see Part A4 table above).
3. Confirm the check script's tolerance is appropriate for floating-point
   double-precision arithmetic (~1e-6 to 1e-4 relative error, never
   tighter).

### C2. Sanity-Check Resolution-Independence Bugs by Running at Two Different Grid Sizes

If a computed flux/quantity looks suspicious, rerun the same RegTest with
a different amr.n_cell z-component. If the result is unchanged, the bug
is almost certainly a hardcoded/stub value (see Part A3).

### C3. Verify Build/Wiring With Grep Before Claiming a Phase Complete

Per Contract R2, every phase's PR description must include actual
grep/diff output confirming:

```
grep -rn "radChoice" Source/DataStructs/ERF_DataStruct.H
grep -rn "ERF_<NewFile>.cpp" CMake/BuildERFExe.cmake Source/Radiation/Make.package
grep -rn "RadiationDiagnostics::append" Source/Radiation/*.cpp
```

### C4. Bit-for-Bit Reproducibility Test

Run the same inputs file with erf.radiation_type = "None" (default)
before and after any new phase's changes; plotfiles must be bit-identical
(no physics regression when radiation is disabled). Mirrors UCM F3.

---

## Quick Checklist for New Radiation Phase

**Before opening PR:**

- [ ] File header with @file, @brief, References section citing relevant literature
- [ ] New *Choice struct members (if any) added to SolverChoice with include and init_params(...) call site — confirmed via grep, snippet pasted in PR description
- [ ] Make.package AND CMake/BuildERFExe.cmake both updated — confirmed via grep, snippet pasted in PR description
- [ ] All new debug prints follow [RAD][PhaseN][ClassOrFile::function] tagged format, IO-rank-only, gated on erf.radiation.v >= 1
- [ ] RADIATION_DIAG: lines always print regardless of verbosity
- [ ] No hardcoded vertical-layer counts or loop bounds — always derived from geom[lev].Domain()
- [ ] Per-column vertical sweeps implemented as ONE kernel per column, not one launch per level
- [ ] No host-side amrex::Print() inside device kernels
- [ ] All inline helper functions defined in the header (never split declaration/definition)
- [ ] New RegTest inputs/sounding files follow Exec/CanonicalTests/ABL/* conventions
- [ ] New RegTest inputs files mirror the full structure of the canonical neutral_abl reference
- [ ] Every new RegTest parameter sanity-checked by hand against known physical ranges
- [ ] Check-script tolerances are floating-point-safe (~1e-6 to 1e-4)
- [ ] RadType::None remains default; all new options gated behind TwoStream check
- [ ] Existing Phase N-1 RegTests re-verified to still pass unchanged
- [ ] No new compiler warnings

---

## Known Issues & Workarounds

### Phase 1 — Bug: SolverChoice/build wiring never connected (Fixed via Phase 1b)

**Issue:** ERF_RadStruct.H and Source/Radiation/* were fully correct
in isolation but never wired into SolverChoice, CMake/BuildERFExe.cmake,
or a real call site — compiled into nothing, diagnostics CSV never
produced.

**Workaround/Fix:** See RAD_DEVELOPMENT.md Phase 1b. Added include +
RadChoice radChoice member + init_params(...) call in SolverChoice;
registered ERF_RadiationDiagnostics.cpp in both build systems; added
ERF_AdvanceTwoStreamRadiation.cpp as the real call site.

**Prevention:** Contract R2 — every future PR must paste grep/diff
confirmation of all three wiring points.

### Phase 1 — Bug: SW cumulative optical depth hardcoded to single layer (Fixed via Phase 1c)

**Issue:** tau_cum = rad_choice.tau_per_layer; in
ERF_AdvanceTwoStreamRadiation.cpp ignored the actual vertical grid
resolution entirely — confirmed because the wrong SW_surface value was
identical across two different domain resolutions (20 vs 64 layers).

**Workaround/Fix:**

```
int n_layers = geom[lev].Domain().length(2);
amrex::Real tau_cum = rad_choice.tau_per_layer * static_cast<amrex::Real>(n_layers);
```

**Prevention:** Contract R3, Part A3/B3 above. Always test a suspicious
flux value at two different grid resolutions before hypothesizing about
missing physics.

### Phase 1 (manual, not coding agent) — Sounding file location/format/realism issues

**Issue:** Shared sounding file at wrong location, wrong column format,
unrealistic tau_per_layer/domain combination, and an overly tight
check-script tolerance were all found and fixed manually by the repo
owner after the coding agent's Phase 1/1b/1c work.

**Workaround/Fix:** See RAD_DEVELOPMENT.md "Manual post-Phase-1
corrections" section for full details.

**Prevention:** Contracts R4, R5. Part A4/A5 checklist items above must
be followed by the coding agent for every future phase to prevent
recurrence.

---

## References

- Source/Radiation/RAD_DEVELOPMENT.md — Phase roadmap and history
- Source/UrbanCanopy/UCM_DEVELOPMENT.md / UCM_MPI_SKILLS.md — structural
  template and cross-module lessons this file mirrors
- Exec/CanonicalTests/ABL/MRF_YSUNew_Enhancements/canonical/{neutral_abl,
  sounding_neutral_abl} — canonical inputs/sounding reference template
- Exec/CanonicalTests/Bomex/inputs_bomex — minimal moisture-enabling
  reference pattern
```

