# Phase 17: Simplified SEB — MultiFab Infrastructure + Noah-MP Passthrough

## Objective

Validate the Phase 17 SEB infrastructure implementation:
- **SEB MultiFabs allocated and populated** from either Noah-MP/LSM passthrough or scalar fallback defaults
- **No prognostic physics yet** — allocation and passthrough only; phase 18 adds diagnostic residual computation
- **Backward compatibility maintained** when `seb_enable=false` (default)

## Test Design

### Baseline Scenario (Default)
- **SEB disabled** (`seb_enable=false`)
- **Expected behavior**: Identical to Phase 16 baseline output
- **Validates**: No regression and full backward compatibility

### Feature-On Scenario
- **SEB enabled** (`seb_enable=true`)
- **No LSM active** — uses scalar fallback defaults
- **Expected behavior**: SEB MultiFabs populated with constant fallback values (deterministic, finite, non-crashing)
- **Validates**: Infrastructure wired correctly with safe fallback path

## Input Files

### `input_sounding` (Atmospheric Profile)
**Purpose**: Initialize atmospheric fields (pressure, temperature, water vapor) on the simulation grid.

**Format** (per line):
```
pressure(Pa)  temperature(K)  water_vapor_mixing_ratio(kg/kg)  [additional fields]
```

**Physical Constraints**:
- **Water vapor mixing ratio (qv)**: Typically 0.003–0.015 kg/kg (3–15 g/kg) for mid-latitude troposphere
  - Near surface: ~8–10 g/kg (initialized at 0.008 kg/kg in this test)
  - Upper troposphere: ~4–6 g/kg (decays with height due to decreasing temperature and moisture holding capacity)
  - Must always be ≥ 0.0 (physically represents mass of water vapor per unit mass of dry air)
  - **Sanity check**: qv should monotonically decrease or remain constant with height; never increase sharply
  
**This Test's Profile**:
- 5 levels (header + 4 data levels)
- Surface (0 m) to 1550 m altitude
- qv profile: 0.008 → 0.008 → 0.006 → 0.004 kg/kg (realistic decay with height)
- Adapted from Phase 12 moist cloud test baseline

### Inputs Files

#### `inputs_seb_disabled`
- **seb_enable = false**
- Baseline case for backward-compatibility regression testing
- Expected output identical to Phase 16 (no new SEB diagnostics)

#### `inputs_seb_enabled`
- **seb_enable = true**
- Tests SEB infrastructure with scalar fallback defaults
- Fluxes populated from constant values (all finite, safe no-op)
- No LSM active; all fields derive from `seb_*_default` parameters

## Validation Criteria

1. **Baseline (SEB disabled)**
   - All radiation diagnostics must be finite
   - Output must match Phase 16 regression baseline exactly (bitwise compatibility)

2. **Feature-on (SEB enabled, no LSM)**
   - All SEB MultiFabs must be allocated and contain finite values
   - Fallback scalar defaults must be properly propagated to every surface column
   - qv field must be finite and physically reasonable throughout initialization and simulation

3. **General**
   - No NaN/Inf in any diagnostic output
   - Heating rates must be nontrivial (non-zero) and decay appropriately with height

## Running the Tests

### Baseline Mode
```bash
erf inputs_seb_disabled
# Check: radiation_diag.dat has finite values, matches Phase 16 baseline
```

### Feature-On Mode
```bash
erf inputs_seb_enabled
# Check: SEB MultiFabs initialized and finite
# Check: Scalar fallback defaults correctly populated
# Check: qv profile reasonable (decreases with height, 0.008→0.004 kg/kg)
```

## Prior Bugfix Note (Part A)

**Issue**: Previous version of `input_sounding` contained unrealistic moisture values (qv=0.0 uniformly with height), which does not represent a physical mid-latitude atmosphere.

**Fix**: Replaced with a vertically varying profile based on Phase 12 moist cloud test sounding:
- **Surface qv**: 0.008 kg/kg (~8 g/kg, typical for moderate humidity)
- **Upper levels**: Decay to 0.004 kg/kg at 1550 m (consistent with exponential moisture decay)
- **Source**: Adapted from `TwoStream_DynamicTau_MoistCloud/input_sounding_phase12_moist` to match this test's 5-level vertical grid

This ensures the initialized atmospheric state is physically meaningful and allows proper validation of moisture-dependent optical depth, cloud fraction, and future SEB residual diagnostics (Phase 18+).

## References

- `Source/Radiation/RAD_DEVELOPMENT.md` — Phase 17 Implementation section
- `Source/DataStructs/ERF_RadStruct.H` — RadChoice SEB parameters
- `Source/ERF.H` — SEB MultiFab vector declarations
