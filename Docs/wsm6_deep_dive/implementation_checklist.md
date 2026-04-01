# WSM6 Implementation Checklist (Derived from Deep Dive)

## Phase 1: Fortran-in-place WSM6 runnable in ERF
- [x] Add `Source/Microphysics/WSM6/` with ERF wrapper class and Fortran interface header.
- [x] Add build wiring in `Exec/Make.ERF` and CMake for WSM6 sources.
- [x] Add `MoistureType::WSM6` dispatch in ERF microphysics selection.
- [x] Implement ERF->WSM6 argument bridge for required state (`qv,qc,qi,qr,qs,qg`, thermodynamics, geometry).
- [x] Implement WSM6 init call and runtime call lifecycle.
- [x] Add accumulation mapping (`rain/snow/graupel`) analogous to Morrison diagnostics.
- [ ] Validate single-step run with `moisture_model=WSM6` (no NaN/abort, non-negative hydrometeors).

## Phase 2: Head-comparable behavior and diagnostics
- [ ] Add optional effective radius diagnostics and radiation coupling arguments.
- [ ] Align boundary/real-width behavior with current head microphysics path.
- [ ] Add regression tests mirroring existing Morrison harness patterns.
- [ ] Add documentation for runtime controls and known differences from WRF defaults.
