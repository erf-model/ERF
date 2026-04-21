# Current status, compiled groups and compile status, last commit hash
- Branch/context: `wsm6-moisture-stubs` in `/home/jmsexton/codes/ERF`.
- Latest local commit hash: `1db6cbe9` (`WSM6 cpp scaffold: complete loop body and drop final ignore_unused (compiles clean)`).
- Build context: `build_no_netcdf_wsm6` configured with `ERF_ENABLE_NETCDF=OFF` and `ERF_ENABLE_WSM6_FORT=ON`.
- Full WSM6-enabled executable build completed in that config (`erf_exec` exists per the phase-1 compare report).
- Runtime status in the same config:
  - WSM6 short run: fails on first advance with `Erroneous arithmetic operation` inside `__mp_wsm6_MOD_mp_wsm6_run`.
  - Morrison short run: passes through step 120.
- Targeted compile groups already verified clean:
  - C++ bridge objects: `ERF_InitWSM6.cpp.o`, `ERF_AdvanceWSM6.cpp.o`, `ERF_UpdateWSM6.cpp.o`
  - Fortran bridge objects: `ERF_module_libmassv.F90.o`, `ERF_mp_radar.F90.o`, `ERF_module_mp_wsm6.F90.o`, `ERF_module_mp_wsm6_isohelper.F90.o`
- The phase-1 compare report records the WSM6 failure as a runtime exception, not a build failure.

# Flag mechanism problem and fix direction
- Problem statement: the useful diagnostic state is split between compile-time selection (`ERF_ENABLE_WSM6_FORT` / `ERF_USE_WSM6_FORT`) and runtime debug output (`erf.microphysics_debug`), so the sidecar must not assume one flag covers both roles.
- Current code-search result: there is no literal `microphysics_debug` symbol in `Source/` or `Exec/`; the flag behavior is described in the handoff notes and likely implemented through runtime parameter plumbing, not as a standalone C++ symbol.
- Fix direction gathered so far:
  - Keep the compile-time bridge toggle separate from the runtime diagnostic level.
  - Gate diagnostics by level, not by ad hoc prints.
  - Use the debug level only to enable controlled host-side tracing and small-case dumps.
  - Preserve the bridge ABI and avoid introducing new behavior into the device-side kernel path while instrumenting.

# Morrison flag pattern findings (what is known; if exact lines unknown, explicitly mark as pending)
- Known from the handoff notes:
  - `erf.microphysics_debug` levels are `0/1/2`.
  - Stable tags used in output are `MPDBG PRE`, `MPDBG POST`, and `MPDBG STEP`.
  - Debug level `2` enables full per-cell dumps for small-test usage.
- Known from the current repo search:
  - No direct `microphysics_debug` symbol was found in `Source/` or `Exec/` with `rg`.
  - This means the exact implementation lines are still pending identification in the current tree.
- Pending:
  - Exact Morrison source lines for the debug flag parsing and the print gates.
  - Exact line numbers for `MPDBG` emission in Morrison and WSM6 source files.
  - Whether the final implementation is purely runtime-param driven or also uses a compile-time guard around the dump path.

# GPU-safe diagnostic print plan details gathered
- Plan already captured in the handoff notes:
  - Use fixed-format `printf` output.
  - Emit `MPDBG PRE`, `MPDBG POST`, `MPDBG STEP` tags so compare output is grep-able.
  - Restrict full per-cell dumps to debug level `2` and small tests.
- GPU-safety direction inferred from the current WSM6 C++ bridge shape:
  - Keep prints on host-side control flow around the kernel, not inside the device lambda.
  - Use guards so the kernel path remains side-effect free.
  - Prefer summary prints first, then full dumps only when the case is tiny and intentionally diagnostic.
- The current WSM6 path already uses host-managed arrays and per-column work buffers in `ERF_AdvanceWSM6.cpp`, which makes host-gated diagnostics the safer insertion point.

# fcompare harness locations and scripts
- Main WSM6 compare workflow README:
  - `Source/Microphysics/WSM6/README`
  - Uses `${INSTALL_DIR}/bin/amrex_fcompare` by default and compares `plt_wsm6_run00120` against `plt_morr_run00120`.
- Stepwise audit script:
  - `Exec/CanonicalTests/SquallLine_2D/run_r600b_subtaskA.sh`
  - Uses `/home/jmsexton/codes/ERF/tools/fcompare_serial`
  - Writes `audit_r600b/fcompare_r600b_steps600_700.tsv`, `audit_r600b/onset_summary.txt`, and raw compare logs.
- Additional WSM6/SquallLine harness artifacts already present:
  - `Exec/CanonicalTests/SquallLine_2D/scripts/run_wsm6_long.sh`
  - `Exec/CanonicalTests/SquallLine_2D/inputs_moisture_WSM6`
  - `Exec/CanonicalTests/SquallLine_2D/inputs_moisture_Morrison_short`
  - `Exec/CanonicalTests/SquallLine_2D/audit_matrix_fcompare.tsv`
  - `Exec/CanonicalTests/SquallLine_2D/audit_matrix_summary.tsv`
  - `Exec/CanonicalTests/SquallLine_2D/audit_r200b/` and `Exec/CanonicalTests/SquallLine_2D/audit_r600b/`
- Direct compare commands currently documented:
  - `amrex_fcompare plt_wsm6_run00120 plt_morr_run00120`
  - `tools/fcompare_serial <pltA> <pltB>`

# Remaining build work (GNUmake, CMake, define locations)
- The WSM6 build toggles are already wired in both build systems:
  - GNUmake define path: `Exec/Make.ERF`
  - GNUmake package path: `Source/Microphysics/WSM6/Make.package`
  - CMake option: `CMakeLists.txt` (`ERF_ENABLE_WSM6_FORT`)
  - CMake source inclusion: `CMake/BuildERFExe.cmake`
  - CMake compiler/option propagation: `CMake/CrayCompilerDetection.cmake`, `CMake/SetAmrexOptions.cmake`
- Current define locations known from the tree:
  - `Exec/Make.ERF` sets `-DERF_USE_WSM6_FORT` when `USE_WSM6_FORT` is `TRUE`.
  - `CMake/BuildERFExe.cmake` sets `target_compile_definitions(... ERF_USE_WSM6_FORT)`.
- Remaining build work is therefore mostly verification and parity, not discovery:
  - Confirm any future WSM6 debug/diagnostic flag additions are mirrored in both GNUmake and CMake paths.
  - Keep the runtime debug parameter and compile-time bridge define separate.
  - If the bridge grows, verify the same include/source lists stay synchronized between `Exec/Make.ERF` and `CMake/BuildERFExe.cmake`.

# Open questions / next steps
- What is the exact first bad arithmetic operation inside `mp_wsm6_run` in the WSM6 short run?
- Which exact source line in the Fortran path throws once the build is switched to a debug/symbolized configuration?
- Is the diagnostic flag plumbing entirely runtime-based, or does Morrison also rely on a compile-time switch for some prints?
- Do any missing GNUmake/CMake parity items remain for the diagnostic-print path beyond the current `ERF_USE_WSM6_FORT` wiring?
- Next useful validation step:
  - Build a debug/symbolized WSM6 executable.
  - Re-run the short WSM6 case.
  - Symbolize the first failing frame and pin the exact line.
