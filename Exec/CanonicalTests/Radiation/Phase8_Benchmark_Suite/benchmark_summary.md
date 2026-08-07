# Phase 8 Benchmark Suite Results
**Timestamp:** 2026-08-07T23:37:36.818235
## Summary
- **Total Cases:** 5- **Passed:** 5- **Failed:** 0- **Overall Status:** ✅ PASS
## Case Results
| Case | Name | Status | Errors |
|------|------|--------|--------|
| `lw_isothermal` | LW isothermal baseline | ✅ PASS | None |
| `phase6_timing` | Coupled SW+LW non-isothermal time-integration | ✅ PASS | None |
| `sw_clearsky` | Clear-sky SW baseline | ✅ PASS | None |
| `sw_cloud_layer` | Cloud-layer absorption | ✅ PASS | None |
| `sw_scattering` | Cloud scattering | ✅ PASS | None |

## Detailed Results
### lw_isothermal: LW isothermal baseline [✅ PASS]
**Status:** No errors

**Key Metrics:**
- `row_count`: 2
- `row_count_pass`: True
- `step_multiplicity_pass`: True
- `callsite_mode_pass`: True
- `LW_surface_mean`: 3.844000e+02
- `LW_surface_final`: 3.844000e+02
- `LW_surface_max`: 3.844000e+02
- `LW_surface_min`: 3.844000e+02
- `LW_TOA_mean`: 3.844000e+02
- `LW_TOA_final`: 3.844000e+02
- `LW_TOA_max`: 3.844000e+02
- `LW_TOA_min`: 3.844000e+02
- `heating_rate_max_mean`: 1.000000e-15
- `heating_rate_max_final`: 1.000000e-15
- `heating_rate_max_max`: 1.000000e-15
- `heating_rate_max_min`: 1.000000e-15
- `heating_rate_max_cv`: 0.000000e+00

### phase6_timing: Coupled SW+LW non-isothermal time-integration [✅ PASS]
**Status:** No errors

**Key Metrics:**
- `row_count`: 10
- `row_count_pass`: True
- `step_multiplicity_pass`: True
- `callsite_mode_pass`: True
- `SW_TOA_mean`: 6.805000e+02
- `SW_TOA_final`: 6.805000e+02
- `SW_TOA_max`: 6.805000e+02
- `SW_TOA_min`: 6.805000e+02
- `LW_surface_mean`: 3.844000e+02
- `LW_surface_final`: 3.844000e+02
- `LW_surface_max`: 3.844000e+02
- `LW_surface_min`: 3.844000e+02
- `heating_rate_max_mean`: 2.015000e-02
- `heating_rate_max_final`: 2.020000e-02
- `heating_rate_max_max`: 2.020000e-02
- `heating_rate_max_min`: 2.010000e-02
- `heating_rate_max_cv`: 2.481390e-03

### sw_clearsky: Clear-sky SW baseline [✅ PASS]
**Status:** No errors

**Key Metrics:**
- `row_count`: 2
- `row_count_pass`: True
- `step_multiplicity_pass`: True
- `callsite_mode_pass`: True
- `SW_TOA_mean`: 6.805000e+02
- `SW_TOA_final`: 6.805000e+02
- `SW_TOA_max`: 6.805000e+02
- `SW_TOA_min`: 6.805000e+02
- `SW_surface_mean`: 3.002500e+02
- `SW_surface_final`: 3.002500e+02
- `SW_surface_max`: 3.002500e+02
- `SW_surface_min`: 3.002500e+02
- `heating_rate_max_mean`: 1.500000e-02
- `heating_rate_max_final`: 1.500000e-02
- `heating_rate_max_max`: 1.500000e-02
- `heating_rate_max_min`: 1.500000e-02
- `heating_rate_max_cv`: 0.000000e+00

### sw_cloud_layer: Cloud-layer absorption [✅ PASS]
**Status:** No errors

**Key Metrics:**
- `row_count`: 2
- `row_count_pass`: True
- `step_multiplicity_pass`: True
- `callsite_mode_pass`: True
- `SW_TOA_mean`: 6.805000e+02
- `SW_TOA_final`: 6.805000e+02
- `SW_TOA_max`: 6.805000e+02
- `SW_TOA_min`: 6.805000e+02
- `SW_surface_mean`: 1.501200e+02
- `SW_surface_final`: 1.501200e+02
- `SW_surface_max`: 1.501200e+02
- `SW_surface_min`: 1.501200e+02
- `heating_rate_max_mean`: 2.500000e-02
- `heating_rate_max_final`: 2.500000e-02
- `heating_rate_max_max`: 2.500000e-02
- `heating_rate_max_min`: 2.500000e-02
- `heating_rate_max_cv`: 0.000000e+00

### sw_scattering: Cloud scattering [✅ PASS]
**Status:** No errors

**Key Metrics:**
- `row_count`: 2
- `row_count_pass`: True
- `step_multiplicity_pass`: True
- `callsite_mode_pass`: True
- `SW_TOA_mean`: 6.805000e+02
- `SW_TOA_final`: 6.805000e+02
- `SW_TOA_max`: 6.805000e+02
- `SW_TOA_min`: 6.805000e+02
- `SW_surface_mean`: 1.200800e+02
- `SW_surface_final`: 1.200800e+02
- `SW_surface_max`: 1.200800e+02
- `SW_surface_min`: 1.200800e+02
- `heating_rate_max_mean`: 3.200000e-02
- `heating_rate_max_final`: 3.200000e-02
- `heating_rate_max_max`: 3.200000e-02
- `heating_rate_max_min`: 3.200000e-02
- `heating_rate_max_cv`: 0.000000e+00

## Tolerance Configuration
- SW_TOA relative error: 0.1%
- SW_surface relative error: 1.0%
- LW_net_surface relative error: 1.0%
- Heating CV upper bound: 5.0%
- Row count tolerance: ±2
