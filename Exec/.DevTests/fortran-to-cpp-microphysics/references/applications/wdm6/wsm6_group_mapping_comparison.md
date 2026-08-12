# WDM6 vs WSM6: Group Map and Variable-to-Tag Alignment

## Overview
This document compares the WDM6 group map (`group_map.md`) structure against the WSM6 variable-to-tag mapping and canonical retreat order defined in `wsm6_implementation_notes.md` (section: Variable-to-tag mapping).

## Structural Comparison

### WDM6 Group Map Structure
The WDM6 group map (lines 20-54 of group_map.md) defines **23 logical process groups**:
- **Prelude/Setup**: P0 (lines 592-659)
- **Minor-step init**: G1a (665-669), G1b (671-677), G1c (683-713)
- **Core physics**: G2-G17 (719-2072)

Each group has:
- **Fortran line range**: exact source location
- **Block purpose**: functional description
- **Native C++ status**: implementation maturity level (partial scaffolding, present, simplified heuristics, not yet ported)

### WSM6 Variable-to-Tag Mapping
The WSM6 implementation notes (lines 671-690) define **variable→tag dependencies**:

| Variable | Associated Tags/Groups |
|----------|----------------------|
| `qr` (rain) | NISLFV_R, PRAUT, PRACW, PREVP |
| `qs` (snow) | NISLFV_SG, PSACI, PRACS, PSAUT, PSEVP |
| `qg` (graupel) | NISLFV_SG, PSACI, PRACS, PSAUT, PSEVP |
| `qi` (ice) | VICE, PRACI, PIDEP, PSAUT |
| `qv`/`q` (vapor) | PREVP, PCOND, QSAT, QSAT2 |
| `t` (temperature) | MELT, PHASE, UPDATE, PCOND |
| `fall_r/s/g` (fallspeeds) | FALL |
| `rho`/`den` (density) | SLOPE1, SLOPE2, SLOPE3 |

### WSM6 Canonical Retreat Order
The validation retreat order (lines 682-686 of wsm6_implementation_notes.md):
```
DENFAC → QSAT → SLOPE1 → NISLFV_R → NISLFV_SG → FALL → SLOPE2 →
MELT → VICE → PRECIP → PHASE → SLOPE3 → PRAUT → PRACW → PREVP →
PRACI → PSACI → PRACS → PSEML → PIDEP → PSAUT → PSEVP → UPDATE →
QSAT2 → PCOND
```

## WDM6 Alignment Strategy

### Key Difference: WDM6 Double-Moment Structure
WDM6 introduces **number concentrations** that WSM6 does not have:
- `nc` (cloud number)
- `nn` (CCN/aerosol number)
- `nr` (rain number)

**Impact on group structure**:
- Groups G1b (DENFAC), G1c (QSAT), and core process groups must handle these additional prognostic variables
- The 3D species-indexed arrays in WSM6 (rslope_r, rslope_s, rslope_g) remain structurally similar in WDM6
- Number-concentration rates must be computed wherever mass rates are, extending the process-group scope

### Mapped Group Correspondence

The WDM6 groups should inherit the WSM6 retreatorder where physics is identical:

| WDM6 Group | WSM6 Tag Equivalent | Variables Affected | WDM6 Divergence |
|-----------|------|-----------|---|
| G1b | DENFAC | `denfac`, density factors | Same structure; uses `nc/nn/nr` in subsequent groups |
| G1c | QSAT | `qs(:,1:2)`, `rh(:,1:2)` | Same for mixing ratios; adds number-conc saturation? |
| G2 | RATES_ZERO | all process-rate arrays | **Extended**: must zero `*nc`, `*nn`, `*nr` rate arrays |
| G3 | CLOUD_SETUP | `rslopec{,2,3}`, `xni` | Same structure |
| G4 | SLOPE1 | rain packing, first slope | Same structure |
| G5a–G5c | NISLFV_R, NISLFV_SG | sedimentation, slopes | Same structure (3-species fallout logic) |
| G6 | SLOPE2 | repack, second slope | Same structure |
| G7–G10e | MELT, VICE, PRECIP, PHASE | melting, freezing, phase transitions | **Extended**: adds number-concentration phase changes |
| G11 | SLOPE3 | third slope, recompute `avedia` | Same structure |
| G12 | WORKDIFF | `diffac`, `venfac` | Same structure |
| G13a | PRAUT (warm-rain rates) | `praut`, `nraut`, `pracw`, `nracw`, etc. | **Extended**: must include `nc` source/sink terms |
| G13b–G13g | PSACI, PSACR, PIDEP, PSAUT, PSEVP | ice/snow/graupel accretion, deposition, auto-conversion, evaporation | **Extended**: number rates for each |
| G14 | UPDATE | state update, latent heat | **Extended**: `nc`, `nn`, `nr` bounds, conservation checks |
| G15 | QSAT2 | recompute `qs(:,1:2)`, `rh(:,1)` | Same structure |
| G16a–G16b | (rain evaporation cleanup, activation) | slope repack, `avedia`, `pcond`, `ncact`, `pcact` | **Critical**: `ncact` (CCN activation) links `nn → nc` |
| G17 | (padding/bounds) | zero tiny condensate, enforce bounds | Same structure |

## Pre/Post Print Requirements by Group

### Core State Variables (All Groups)
**Pre-group inputs** (computed in prior groups or supplied):
- Density: `den`, `den_tmp` (from P0 or sedimentation)
- Temperature: `t`, thermodynamic coefficients `xl`, `cpm`
- Mixing ratios: `q`, `qc`, `qi`, `qr`, `qs`, `qg`
- **WDM6 addition**: `qnc`, `qnn`, `qnr` (or `nc`, `nn`, `nr`)
- Intermediate slopes/saturation: `rslope_r/s/g`, `qsat_r/s/g`, `rh_r/s/g`
- Sedimentation parameters: `fall_r/s/g`, `work1_r/s/g`

**Post-group outputs** (what changed in the group):
- Updated mixing ratios and number concentrations
- Updated process-rate arrays (mass and number)
- Updated sedimentation/fallout accumulations

### Group-Specific Pre/Post Variables

#### G1b (DENFAC)
| Requirement | Variable | Notes |
|-------------|----------|-------|
| **PRE** | `vrec_d`, `vsqrt_d` | density-dependent factors |
| **POST** | `denfac` | output for downstream groups |

#### G1c (QSAT)
| Requirement | Variable | Notes |
|-------------|----------|-------|
| **PRE** | `t`, `p`, `den` | temperature, pressure, density |
| **POST** | `qs(:,1:2)`, `rh(:,1:2)` | saturation and RH for liquid/ice |
| **WDM6 extension** | `nn`, `ccn0`, `xland` | affects activation pathway |

#### G2 (RATES_ZERO)
| Requirement | Variable | Notes |
|-------------|----------|-------|
| **PRE** | none (zeroing operation) | — |
| **POST** | all `p*` arrays + `n*` arrays (WDM6) | praut, pracw, pigen, pidep, psaut, pgaut, nraut, ncact, nnact (WDM6), etc. |

#### G13a (PRAUT - Warm-Rain Rates)
| Requirement | Variable | Notes |
|-------------|----------|-------|
| **PRE** | `q`, `qc`, `t`, `den`, `denfac`, `qs`, `rh` | state and derived thermo |
| **PRE_WDM6** | `nc`, `nn` | cloud and CCN number for nucleation |
| **POST** | `praut`, `pracw`, `prevp` (mass), `nraut`, `nracw`, `nccol`, `nrcol` (number) | rain auto, cloud accretion, evap + their number rates |
| **POST_WDM6** | `ncact`, `pcact` (CCN activation), `nn` sink | WDM6-specific CCN activation pathway |

#### G14 (UPDATE)
| Requirement | Variable | Notes |
|-------------|----------|-------|
| **PRE** | all `p*` rates (mass), all `n*` rates (number, WDM6) | integrated process changes |
| **PRE** | `xl`, `cpm`, `q`, `qc`, `qi`, `qr`, `qs`, `qg` | old state |
| **PRE_WDM6** | `nc`, `nn`, `nr` | old number concentrations |
| **POST** | updated `t`, `q`, `qc`, `qi`, `qr`, `qs`, `qg` | new state after rate application + limiting |
| **POST_WDM6** | updated `nc`, `nn`, `nr` + bounds checks | new numbers + conservation enforcement |

#### G16b (PCOND + Activation)
| Requirement | Variable | Notes |
|-------------|----------|-------|
| **PRE** | `q`, `qc`, `t`, `qs`, `rh`, all prior process outputs | condensation/evaporation inputs |
| **PRE_WDM6** | `nn`, `ccn0`, `xland` | aerosol activation inputs |
| **POST** | `pcond`, `ncact`, `pcact` (activation mass/number), final `qc`, `q` | condensation + WDM6 CCN-to-cloud link |
| **POST_WDM6** | `nc` updated via activation (`ncact`), `nn` consumed | CCN → cloud number feedback |

## Validation Strategy Implications

### 1. WSM6 Retreat Order Applicability
The WSM6 canonical retreat order can be **directly used for WDM6 dry/liquid-only regimes**:
```
DENFAC → QSAT → SLOPE1 → NISLFV_R → NISLFV_SG → ...
```

### 2. WDM6 Extension Points (Must be Added)
- **G2 number-array zeroing**: extend RATES_ZERO to include `*nc`, `*nn`, `*nr` arrays
- **G13a activation logic**: add `ncact`, `pcact`, `nnact` (WDM6-specific CCN activation)
- **G14 number updating**: add bounds-checking and conservation for number concentrations
- **G16b activation feedback**: link `nn` consumption and `nc` source via `ncact`

### 3. Pre/Post Print Tags for WDM6 Validation
When running WDM6 in bridge or native mode, print before/after each group:

**Minimal tag set (fast retreat)**:
- Before G2: `PRE_RATES_ZERO` (state snapshot)
- After G2: `POST_RATES_ZERO` (zeroed arrays)
- Before G13a: `PRE_PRAUT` (`qc`, `denfac`, `nc`, `nn`, `qs`, `rh`)
- After G13a: `POST_PRAUT` (`praut`, `pracw`, `nraut`, `ncact`, `nnact`)
- Before G14: `PRE_UPDATE` (all process rates, old state)
- After G14: `POST_UPDATE` (new state, bounds-checked numbers)

**Extended tag set (narrow retreat)**:
- Add line-level forensic prints inside multi-rate groups (G13a, G13b, G14)
- Use per-variable `PRE_<function>` / `POST_<function>` tags around key activation/phase-change function calls

## Known WDM6-Specific Divergences

From `wdm6.md` and `group_map.md`:

1. **CCN/aerosol handling**: `nn` (CCN number) must be initialized from `wdm6.ccn0` and preserved through state-copy; WSM6 has no analogous variable
2. **Cloud number activation**: G16b must link `nn` → `nc` via `ncact`/`pcact`; WSM6 computes `ncact` but lacks the aerosol pool
3. **Xland dependency**: `xland` (land fraction) affects activation pathway; not present in WSM6
4. **Missing ice double-moment**: WDM6 is liquid-only double-moment; future variants may add ice number

## Recommended Next Steps

1. **Extract WDM6-specific pre/post variables for each group** into a TSV parallel to WSM6's validation_manifest structure
2. **Implement G2 extension** (number array zeroing) as first bounded slice validation
3. **Implement G13a extension** (warm-rain rates including CCN activation) as early strong frontier candidate
4. **Establish WDM6 Bubble campaign** with full pre/post tag logging before formal validation push
5. **Link validation_manifest rows** to this comparison table so operators know which variables to monitor per group

