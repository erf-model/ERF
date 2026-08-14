# WDM6 Group Map

Canonical file-order group map for the frozen tracked Fortran source:
`Source/Microphysics/WDM6/ERF_module_mp_wdm6.F90`.

This is the working WDM6 analogue of the WSM6 group/frontier structure.
Line numbers below refer to the tracked in-repo file, not the untracked
top-level `module_mp_wdm6.F` comment copy. The top-level file is guidance
for naming and block intent only.

## Scope

- Driver pre/post wrapper: `wdm6`
- Physics kernel to port in bounded slices: `wdm62D`
- Active kernel span: lines 560-2075 in
  `ERF_module_mp_wdm6.F90`

## File-Order Groups

The first four columns describe the frozen Fortran and are current by
construction: they change only when the tracked source does.

The fifth column does **not** track current state and must not be read as though
it does. It is a first-port impression, written when each group was first
surveyed and never systematically refreshed. Port status is owned by
`campaign_decisions.tsv`, per the ledger model in `validation-operators.md`
(`validation_manifest.tsv` for frontier state, `decisions.tsv` for the mutations
that produced it). Two sources of truth for one fact is a drift generator, and
this column lost the race: by 14 August 2026 it was stale on five rows whose
defects had been closed, and on one row it contradicted the manifest outright.

To answer "where does this group actually stand", grep the group id in
`campaign_decisions.tsv` and read the latest `after_status`. Do not refresh the
column below to match; leave it as the historical impression it is.

| Group | Fortran lines | Block | Purpose | First-port impression (NOT current — see `campaign_decisions.tsv`) |
| --- | --- | --- | --- | --- |
| `P0` | 592-659 | Prelude/setup | Clamp incoming state, derive `dend`, `cpm`, `xl`, `qcr`, `delz_tmp`, `den_tmp`, reset precip accumulators, compute `loops/dtcld` | Partial scaffolding present; not grouped explicitly |
| `G1a` | 665-669 | Minor-step init | Reset `mstep`, `mnstep`, `flgcld` for the active minor loop | Present implicitly |
| `G1b` | 671-677 | `DENFAC` | Compute `denfac` from `vrec_d`/`vsqrt_d` | Present and close to exact |
| `G1c` | 683-713 | `QSAT` | Compute `qs(:,1:2)` and `rh(:,1:2)` | Present and close to exact |
| `G2` | 719-775 | `RATES_ZERO` | Zero all process-rate, fallout, and number-rate work arrays | Present, narrower than Fortran surface |
| `G3` | 777-795 | `CLOUD_SETUP` | Compute `rslopec{,2,3}` and `xni` | Only partial/provisional in native path |
| `G4` | 800-809 | `SLOPE1` | Pack `qrs_tmp/ncr_tmp` and call first `slope_wdm6` | Not yet ported faithfully |
| `G5a` | 813-825 | Rain sedimentation setup | Determine rain/number substep count from `work1(:, :, 1)` and `workn` | Not yet ported faithfully |
| `G5b` | 827-872 | `NISLFV_R` | Rain + rain-number sedimentation with repeated `slope_rain` refresh | Not yet ported faithfully |
| `G5c` | 876-902 | `NISLFV_SG` | Snow/graupel sedimentation via `nislfv_rain_plm6` and lower-boundary slab fall | Not yet ported faithfully |
| `G6` | 903-912 | `SLOPE2` | Repack and call second `slope_wdm6` after sedimentation | Not yet ported faithfully |
| `G7` | 914-969 | `MELT` | Warm-phase snow/graupel melting (`psmlt`, `pgmlt`) plus rain-number feedback | Present only as simplified heuristics |
| `G8` | 973-1001 | `VICE` | Ice fall speed + `nislfv_rain_plmr` ice sedimentation | Present only as simplified heuristics |
| `G9` | 1006-1035 | `PRECIP` | Surface precipitation accumulation and `sr` update | Present only as simplified surface accumulation |
| `G10a` | 1041-1055 | `PHASE` | `pimlt` / `nimlt` cloud-ice melt to cloud water | Present as simplified heuristics |
| `G10b` | 1060-1069 | `PHASE` | `pihmf` homogeneous freezing of cloud water | Present as simplified heuristics |
| `G10c` | 1074-1091 | `PHASE` | `pihtf` heterogeneous cloud-water freezing | Present as simplified heuristics |
| `G10d` | 1096-1113 | `PHASE` | `pgfrz` rain freezing to graupel | Present as simplified heuristics |
| `G10e` | 1117-1122 | Phase cleanup | Clamp `ncr(:, :, 2:3)` non-negative after phase changes | Present implicitly |
| `G11` | 1127-1159 | `SLOPE3` | Repack, third `slope_wdm6`, recompute `avedia` and `rslopec` | Not yet ported faithfully |
| `G12` | 1162-1167 | `WORKDIFF` | Compute `diffac` and `venfac` work arrays | Not yet ported faithfully |
| `G13a` | 1178-1277 | Warm-rain rates | `praut`, `nraut`, `pracw`, `nracw`, `nccol`, `nrcol`, `prevp` | Present but still heuristic |
| `G13b` | 1292-1376 | Ice/rain and ice/snow/graupel accretion | `praci`, `piacr`, `niacr`, `psaci`, `pgaci` | Not yet ported faithfully |
| `G13c` | 1381-1428 | Cloud-water accretion by frozen precip | `psacw`, `nsacw`, `pgacw`, `ngacw`, `paacw`, `naacw` | Not yet ported faithfully |
| `G13d` | 1433-1537 | Rain/snow/graupel collection and enhanced melting | `pracs`, `psacr`, `nsacr`, `pgacr`, `ngacr`, `pgacs`, `pseml`, `nseml`, `pgeml`, `ngeml` | Not yet ported faithfully |
| `G13e` | 1543-1599 | Ice deposition/nucleation | `pidep`, `psdep`, `pgdep`, `pigen` with `ifsat` gating | Not yet ported faithfully |
| `G13f` | 1605-1617 | Cold autoconversion | `psaut`, `pgaut` | Not yet ported faithfully |
| `G13g` | 1624-1641 | Warm evaporation of melting snow/graupel | `psevp`, `pgevp` | Not yet ported faithfully |
| `G14` | 1650-1899 | `UPDATE` | Mass limiting, number limiting, state update, latent-heat update; split by `t<=t0c` and `t>t0c` | Not yet ported faithfully |
| `G15` | 1904-1932 | `QSAT2` | Recompute `qs(:,1:2)` and `rh(:,1)` after `G14` | Not yet ported faithfully |
| `G16a` | 1934-1967 | Rain evaporation cleanup | Repack slope, compute `avedia(:,2)`, apply `di82` rain-to-cloud/CCN collapse | Not yet ported faithfully |
| `G16b` | 1969-2040 | `PCOND` + activation | `ncact`, `pcact`, final `pcond`, cloud evaporation returning `nc -> nn` | Not yet ported faithfully |
| `G17` | 2045-2072 | Padding/bounds | Zero tiny condensate; enforce rain/cloud slope bounds by adjusting `ncr` | Not yet ported faithfully |

## Bridge-wrapper sites (outside the Fortran kernel span)

These are **not** file-order groups and must never be added to the `P0..G17`
sequence: they have no counterpart in `ERF_module_mp_wdm6.F90`, so the
`Fortran lines` column is `NA` by definition. They live in the ERF bridge
caller, outside the `wdm6` driver wrapper and outside the 560-2075 kernel span,
which is exactly why no `PRE_G*`/`POST_G*` pair can observe them. They are the
sites `bubble_campaign.md` directs you to check **before** group retreat.

Numbering is declaration order, not call order.

| Site | Fortran lines | Location | Purpose | Instrumented |
| --- | --- | --- | --- | --- |
| `W1_THETAWB` | NA | `ERF_AdvanceWDM6.cpp`, after `mp_wdm6_run_c` | T -> theta writeback (`theta = t / exner`). Fortran branch only; the native branch never writes `mic_fab_vars[theta]` | yes |
| `W2_PACKUNPACK` | NA | `ERF_module_mp_wdm6_isohelper.F90`, `mp_wdm6_run_c` | Pre-kernel pack / post-kernel unpack across the storage/active bounds split | no |
| `W3_STATECOPY` | NA | `ERF_UpdateWDM6.cpp` | `Copy_Micro_to_State` / `Copy_State_to_Micro` | no |
| `W4_LEGACY` | NA | `ERF_AdvanceWDM6.cpp`, native branch | Retained `#if 0` pre-group regions: Steps 4-8, 9, 10, 11 | n/a, compiled out |

Tag surface, same shape and tier as a canonical group tag but in its own
namespace so it cannot collide with `P0..G17` or be mis-ingested by
`tools/build_ci` as a group closure:

```
WDM6-<FORT|CPP>_<PRE|POST>_<Wn>_<NAME> k <field...>
```

Emit the full column `kts..kte` at `erf.micro_diag_target_column`, per the WSM6
all-k print contract. Gate at `erf.microphysics_debug >= 1`, the Tier 1 gate:
these are forward-validation surfaces, not line-level forensics. Escalate to
Tier 2 only to narrow *within* a site once its pair shows a divergence, and note
that WDM6 has no `micro_diag_mode` plumbing today, so the `WSM6-DIAG-T2`
forensic gate does not yet exist here.

`W1_THETAWB` field order is contractual: `k theta t p exner t_over_exner`.
`t_over_exner` is the decisive column — it is what `theta` should equal if the
writeback ran.

## Port Order

Use these groups in strict file order for bounded native slices:

`P0 -> G1a -> G1b -> G1c -> G2 -> G3 -> G4 -> G5a -> G5b -> G5c -> G6 -> G7 -> G8 -> G9 -> G10a -> G10b -> G10c -> G10d -> G10e -> G11 -> G12 -> G13a -> G13b -> G13c -> G13d -> G13e -> G13f -> G13g -> G14 -> G15 -> G16a -> G16b -> G17`

## Immediate Use

For the next bounded implementation slice, prefer the earliest group that is:

- still not faithful in the native path
- structurally upstream of most later divergence
- bridge-verifiable with a compact tag/block signature

At the time this map was written, the first practical candidate is
`G3` or `G13a`, depending on whether the validation surface already
exposes enough pre-`praut` state to separate setup errors from warm-rain
rate errors.
