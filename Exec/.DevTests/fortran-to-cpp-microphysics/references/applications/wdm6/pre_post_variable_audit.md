# WDM6 PRE/POST Diagnostic Variable Audit (Groups G3-G11)

Audit of `Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp` (clean baseline, post-stash of
forensic G10E/G11 WIP prints) against the required variable sets described in
`wsm6_group_mapping_comparison.md`. Scope: implemented groups G3 through G11 (the
groups that currently carry `WDM6-CPP_PRE_G*`/`WDM6-CPP_POST_G*` tags). G1b/G1c/G2/G12+
are intentionally out of scope here since they have no tags yet in the native path.

Line numbers below refer to the clean baseline file (2264 lines total).

## G3 — CLOUD_SETUP (~927-971)
- PRE: `rslopec, rslopec2, rslopec3, xni, qc, nc, qi, den` (8) — matches POST exactly.
- POST: same 8 values.
- **Status: OK.** No gaps.

## G4 — SLOPE1 (~977-1055)
- PRE: `qr, qs, qg, nr, denfac, t, den` (7).
- POST: `rslope[0,1,2], rslopeb[0,1,2], work1[0,1,2], workn, qr, qs, qg, nr, denfac, t, den` (17).
- **Gap:** `rslope2[0..2]` and `rslope3[0..2]` (squared/cubed slope terms) are computed by
  `slope_wdm6` but not printed in POST_G4.

## G5a — Rain sedimentation setup (~1066-1115)
- PRE: `work1(...,0), workn, delz, dtcld` (4).
- POST: `work1(...,0), workn, mstep, mstepmax, numdt, delz, dtcld` (7).
- **Status: OK** — POST is a superset of PRE plus the derived substep counters; no known gaps.

## G5b — NISLFV_R (~1122-1227)
- PRE: `qr, nr, work1(...,0), workn, mstep, delz` (6).
- POST: `qr, nr, rslope[0], rslopeb[0], rslope2[0], rslope3[0], work1[0], workn, delz, dtcld` (10).
- **Status: OK.** Rain-only slope terms (index 0) are all present in POST.

## G5c — NISLFV_SG (~1233-1306)
- PRE: `qs, qg, worka(inline), den*qs, den*qg, delz, dtcld` (7).
- POST: `qs, qg, worka(recomputed), den*qs, den*qg, work1(klo,1), work1(klo,2)` (7).
- **Status: OK.** No known gaps for the sedimentation-fallout values this group owns.

## G6 — SLOPE2 (~1312-1390)
- PRE: `qr, qs, qg, nr` (4).
- POST: `rslope[0,1,2], work1[0,1,2]` (6).
- **Gap:** `rslopeb[0..2]`, `rslope2[0..2]`, `rslope3[0..2]`, and `workn` are all computed by
  the second `slope_wdm6` call but omitted from POST_G6 (same class of gap as G4).

## G7 — MELT (psmlt/pgmlt) (~1393-1488)
- PRE: `qr, qs, qg, nr, t` (5).
- POST: `qr, qs, qg, nr, t` (5) — confirmed via direct read at line 1479; identical variable
  set to PRE, i.e. only final state is printed.
- **Gap:** intermediate melting rates `psmlt`, `pgmlt` (and the rain-number feedback term
  `gfac`) are not printed in either PRE or POST, so the melting-rate contribution itself is
  not independently verifiable from the diagnostic stream — only the net state change is.

## G8 — VICE (ice fall speed + sedimentation) (~1494-1568)
- PRE: `qi, xni, den, t` (4) — confirmed at line 1496.
- POST: `qi, fallc(klo), den*qi` (3) — confirmed at line 1561.
- **Gap:** `xni` and `t` (both PRE inputs) are not echoed in POST, and the intermediate
  `work1c` (ice fall speed) is not printed at all, so the fall-speed computation itself
  cannot be checked independently of the sedimentation-column routine's net effect.

## G9 — PRECIP (~1573-1628)
- PRE: `work1(klo,0)[fall_r], work1(klo,1)[fall_s], work1(klo,2)[fall_g], fallc(klo)[fall_c], dtcld` (5).
- POST: `rain_arr(klo), sr_arr(klo), snow_arr(klo), graup_arr(klo)` (4).
- **Status: OK.** Inputs and accumulator outputs are both fully captured.

## G10a — PHASE (cloud-ice melt) (~1634-1674)
- PRE/POST: `qc, qi, nc, xni, t` (5, identical).
- **Status: OK.** No gaps.

## G10b — PHASE (homogeneous freezing) (~1770-1813)
- PRE/POST: `qc, qi, nc, t, (t0c-t)[supcol]` (5, identical).
- **Status: OK.** No gaps.

## G10c — PHASE (heterogeneous freezing, pihtf) (~1819-1879)
- PRE/POST: `qc, qi, nc, t, (t0c-t)` (5, identical).
- **Gap:** intermediate rates `pfrzdtc`, `nfrzdtc`, and the cubed cloud-droplet slope term
  `rs3` (= `rslopec^3`) are computed but never printed, so the freezing-rate magnitude is
  only inferable from before/after `qc`/`qi`/`nc` deltas.

## G10d — PHASE (rain-to-graupel freezing, pgfrz) (~1885-1942)
- PRE/POST: `qr, qg, nr, t, (t0c-t)` (5, identical).
- **Gap:** intermediate rates `pfrzdtr`, `nfrzdtr`, and `rs3` (= `rslope3_arr(...,0)`) are
  computed but not printed — same class of gap as G10c.

## G10e — Phase cleanup (clamp nc/nr) (~1947-1970)
- PRE/POST: `nc, nr` (2, identical).
- **Status: OK** for what this group does (a pure clamp on two variables); narrow by design,
  not a gap relative to the group's actual scope.

## G11 — SLOPE3 (~1976-2094)
- PRE: `qr, qs, qg, nr` (4).
- POST: `rslope[0][rain only], avedia[1], rslopec, rslopec2, avedia[0], work2` (6).
- Kernel actually computes, beyond what's printed:
  - `rslope[1,2]` (snow/graupel), `rslopeb[0,1,2]`, `rslope2[0,1,2]`, `rslope3[0,1,2]`,
    `work1[0,1,2]` (rain/snow/graup terminal velocities), `workn` — from the repacked
    slope_wdm6 calls (G11b).
  - `rslopec3` — from the lamdac recompute (G11c), alongside `rslopec`/`rslopec2` which
    *are* printed.
  - `work1[0]`/`work1[1]` (diffac results for liquid/ice) — from G11d; only `work2`
    (venfac) is printed, not the two diffac outputs that feed it conceptually.
  - Updated state vars `qr, qs, qg, nr, t, den` are not echoed in POST (only carried
    forward implicitly).
- **Gap:** POST_G11 omits snow/graupel `rslope[1,2]`, `rslopeb[0..2]`, `rslope3[0..2]`,
  `rslopec3`, `work1[0]`/`work1[1]` (diffac outputs), and the state variables. This is the
  largest gap of any implemented group, consistent with the file's own comment marking the
  disabled block after G10a as noting G11 is a "slope repack" step still needing broader
  validation coverage.

## Untagged computational steps (no PRE/POST tags at all)
These correspond to comparison-doc groups G1b/G1c/G2 but currently have zero diagnostic
tags in the native path:
- **"Step 1" (denfac calc, ~line 867)** — no tag. Comparison doc requires PRE=`vrec_d,vsqrt_d`,
  POST=`denfac`.
- **"Step 2" (saturation calc qsatw/qsati/rhw/rhi, ~874-904)** — no tag. Comparison doc
  requires PRE=`t,p,den`, POST=`qs(:,1:2),rh(:,1:2)`, WDM6 extension=`nn,ccn0,xland`.
- **"Step 3" (zero process rates, ~909-919)** — no tag. Comparison doc requires
  POST=all `p*`/`n*` rate arrays.

These are out of scope for "add missing fields to *existing* prints" since no print exists
yet — flagged here for visibility but not touched by the edits below.

## Summary of Edits to Make (existing prints only)
1. **POST_G4**: add `rslope2[0,1,2]`, `rslope3[0,1,2]`.
2. **POST_G6**: add `rslopeb[0,1,2]`, `rslope2[0,1,2]`, `rslope3[0,1,2]`, `workn`.
3. **POST_G10C**: add `rs3` (trivially recomputable from `rslopec_arr`; `pfrzdtc`/`nfrzdtc`
   are lambda-local and not persisted, so left for a future scratch-array change).
4. **POST_G10D**: add `rs3` (recomputable from `rslope3_arr(...,0)`; `pfrzdtr`/`nfrzdtr`
   left for the same reason as G10C).
5. **POST_G11**: add `rslope[1,2]` (snow/graup), `rslopeb[0,1,2]`, `rslope3[0,1,2]`,
   `rslopec3`, `work1[0]`, `work1[1]`.

Not addressed here (would require restructuring rather than adding a field to an existing
printf, or require capturing values not currently stored past their local scope):
- G7 POST melting rates (`psmlt`/`pgmlt` are local to the `ParallelFor` lambda body, not
  stored in a persistent array — would need a scratch array to expose them at PRE/POST
  scope outside the kernel).
- G8 POST `work1c` (also lambda-local in its own separate `ParallelFor`, not persisted).
- G10c/G10d `pfrzdtc`/`nfrzdtc`/`pfrzdtr`/`nfrzdtr` are also lambda-local — same
  constraint as G7/G8. `rs3` is trivially recomputable at print time (from
  `rslopec_arr`/`rslope3_arr`), so it is added; the rate terms are not since they'd need a
  scratch FArrayBox to survive past the kernel that computes them.
