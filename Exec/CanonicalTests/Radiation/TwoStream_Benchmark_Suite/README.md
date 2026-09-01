# Phase 8: Validation & Benchmarking Suite for TwoStream Radiation

## Overview

Phase 8 establishes a canonical validation and benchmarking workflow for the TwoStream radiation module. This suite provides:

- **Repeatable test matrix** covering SW, LW, clouds, and coupled cases
- **Automated metric extraction** from diagnostic CSV files
- **Central tolerance configuration** for pass/fail thresholds
- **Machine-readable and human-readable reports** (JSON + Markdown)
- **CI-friendly behavior** with non-zero exit codes on failure

## Benchmark Cases

The suite includes 5 benchmark cases, each testing a specific aspect of the radiation model:

| Case | Name | Physics | Expected Steps | Rows/Step | Mode |
|------|------|---------|----------------|-----------|------|
| `sw_clearsky` | Clear-sky SW baseline | Beer-Lambert direct-beam solar | 1 | 2 | both |
| `lw_isothermal` | LW isothermal baseline | Gray-gas with energy balance check | 1 | 2 | both |
| `sw_cloud_layer` | Cloud-layer absorption | Cloud optical depth attenuation | 1 | 2 | both |
| `sw_scattering` | Cloud scattering | Two-stream scattering approximation | 1 | 2 | both |
| `phase6_timing` | Coupled SW+LW time-stepping | RhoTheta coupling over 10 steps | 10 | 1 | pre_only |

### Case Descriptions

#### Case 1: Clear-sky SW baseline
- **Physics:** Beer-Lambert solar radiation with direct-beam attenuation
- **Config:** `../SW_ClearSky_Analytical/inputs`
- **Parameters:** S₀=1361 W/m², zenith=60°, τ=0.05/layer
- **Key metrics:** SW_TOA should match analytical value, flux positive
- **Expected behavior:** Deterministic flux at TOA and surface

#### Case 2: LW isothermal baseline
- **Physics:** Gray-gas longwave radiation in isothermal mode
- **Config:** `../LW_Isothermal/inputs`
- **Parameters:** T_iso=288.15K, τ_lw=1.0/layer
- **Key metrics:** Heating rate should be zero (energy balance test)
- **Expected behavior:** Net heating is exactly zero in isothermal mode

#### Case 3: Cloud-layer absorption
- **Physics:** Shortwave with cloud layer optical depth
- **Config:** `../SW_Cloud_Layer/inputs`
- **Parameters:** Height-varying cloud τ values
- **Key metrics:** Surface flux reduced by cloud absorption
- **Expected behavior:** Clear sky → cloud layer → attenuated surface

#### Case 4: Cloud scattering
- **Physics:** Meador-Weaver two-stream approximation for scattering
- **Config:** `../SW_Scattering_Cloud/inputs`
- **Parameters:** Cloud layer with single-scattering albedo
- **Key metrics:** Surface flux further reduced by scattering
- **Expected behavior:** Scattering increases back-scattering vs. pure absorption

#### Case 5: Coupled SW+LW time-stepping
- **Physics:** Coupled SW + LW radiation with RhoTheta heating over 10 steps
- **Config:** `../TwoStream_TimeIntegration/inputs` (modified for pre_only mode)
- **Parameters:** dt=0.5s, stop_time=5.0s (10 steps), diag_callsite_mode=pre_only
- **Key metrics:** Stable heating rates, no NaN/Inf, correct row count for pre_only
- **Expected behavior:** Demonstrates time-stepping consistency and diagnostics mode control

## Running the Benchmark Suite

### Prerequisites
- ERF executable (`erf.ex`) available and in `$PATH` or specified
- Radiation module compiled with diagnostics support
- Python 3.6+

### Quick Start

From the `TwoStream_Benchmark_Suite` directory:

```bash
# Run validation on existing outputs (assumes cases have been run)
python3 run_benchmark_suite.py --verbose

# Run with minimal output
python3 run_benchmark_suite.py
```

### Complete Workflow (with test execution)

**Note:** The current version of `run_benchmark_suite.py` validates outputs from existing case directories. To run the full benchmark suite with test execution, you would:

1. Build ERF with radiation support
2. Run each test case individually or via integration with CI/CD

Example for single case (manual):
```bash
cd ../SW_ClearSky_Analytical
mpirun -np 1 erf.ex inputs
cd ../../TwoStream_Benchmark_Suite
python3 run_benchmark_suite.py --verbose
```

### Output Artifacts

After running validation, the following files are generated in `TwoStream_Benchmark_Suite/`:

- **`benchmark_summary.json`** — Machine-readable results with:
  - Overall pass/fail status
  - Per-case status and errors
  - Extracted metrics (flux, heating, cadence)
  - Timestamp

- **`benchmark_summary.md`** — Human-readable report with:
  - Summary table of all cases
  - Detailed results for each case
  - Error explanations
  - Tolerance configuration reference

## Metrics and Validation

### Flux Metrics
For each case, the suite extracts:
- **SW_TOA**: Top-of-atmosphere shortwave flux [W/m²]
- **SW_surface**: Surface shortwave flux [W/m²]
- **LW_surface**: Surface longwave net flux [W/m²]

For each flux metric, computed:
- Mean and final values
- Max/min values
- Finite check (no NaN/Inf)

### Heating Metrics
- **heating_rate_max**: Maximum heating rate [K/s] per diagnostic record
  - Mean, final, max, min
  - Coefficient of variation (CV) for stability check
  - Nonzero check (should be > 1e-12 K/s where physically expected)

### Cadence Metrics
- **Row count**: Total diagnostic records (vs. expected)
- **Rows per step**: Should match diag_callsite_mode expectation
- **Call-site distribution**: Validate pre/post/both mode filtering

## Tolerance Configuration

All tolerances are centrally defined in `benchmark_tolerances.py`. See that file for details. Key tolerances:

| Metric | Tolerance | Unit | Rationale |
|--------|-----------|------|-----------|
| SW_TOA relative error | 0.1% | — | Strict: should match analytical solution |
| SW_surface relative error | 1.0% | — | Moderate: includes numerical precision |
| LW_net_surface relative error | 1.0% | — | Moderate: energy balance check |
| Heating rate CV | 5% | — | Stability: allow ±5% variation across steps |
| Row count deviation | ±2 | rows | Startup/teardown variation |
| Heating rate nonzero threshold | 1e-12 | K/s | Machine precision check |

## Diagnostics Modes

The suite tests multiple diagnostics cadence modes to verify Phase 6/7 functionality:

- **`both`** (Cases 1–4): Both pre_dycore and post_dycore calls logged
  - Expected: 2 rows per timestep
  - Validates: Call-site filtering works for pre/post separation

- **`pre_only`** (Case 5): Only pre_dycore calls logged
  - Expected: 1 row per timestep
  - Validates: Single-mode filtering works correctly

This ensures the diagnostic system respects runtime configuration without breaking backward compatibility.

## Exit Codes

The benchmark suite returns:
- **0**: All cases passed ✅
- **1**: One or more cases failed ❌

This enables CI/CD integration for automated testing.

## Troubleshooting

### Missing Diagnostic File
If a case fails with "Diagnostic CSV not found":
- Verify case was run successfully
- Check that `erf.radiation.diag_file` parameter matches expected filename
- Ensure radiation module is enabled in case configuration

### Row Count Mismatch
If row count doesn't match expected:
- Check case configuration (dt, stop_time)
- Verify diag_callsite_mode matches case definition
- Ensure diagnostics aren't being suppressed (diag_enable=true)

### CV Stability Check Failed
If heating rate coefficient of variation is too high:
- May indicate numerical instability or transient effects
- Compare mean vs. final values for drift
- Check for NaN/Inf values (detected separately)

### Call-site Mode Validation Failed
If call-site filtering isn't working:
- Verify `diag_callsite_mode` parameter in inputs file
- Check that diagnostic code correctly tags pre/post calls
- Ensure Phase 6/7 diagnostics module is active

## Phase 6/7 Integration

The benchmark suite preserves and validates Phase 6/7 diagnostics semantics:

✅ **call_site support**: Each record tagged with pre/post identifier
✅ **mode-aware cadence**: `both`, `pre_only`, `post_only` modes respected
✅ **dedup identity**: Not weaker than `(step,time,call_site)` 
✅ **GPU safety**: No host I/O in device code (unchanged from Phase 6/7)

Case 5 specifically tests single call-site mode to ensure Phase 6/7 filtering works.

## Documentation

- **`RAD_DEVELOPMENT.md`**: Phase 8 section with architecture overview
- **`RAD_MPI_SKILLS.md`**: Lesson on benchmark reproducibility and diagnostics-aware validation
- **`README.md`** (this file): User guide and case descriptions

## Future Extensions (Phase 9+)

Potential enhancements:
- Performance profiling and timing metrics
- GPU vs. CPU comparison matrices
- MPI scaling studies
- Regression detection with historical baseline storage
- Automated report visualization
