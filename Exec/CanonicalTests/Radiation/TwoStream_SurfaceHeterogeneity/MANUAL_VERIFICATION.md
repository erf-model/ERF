# Phase 11 Surface Heterogeneity Implementation - Manual Verification Guide

This document provides step-by-step instructions for validating the Phase 11 implementation on your local system.

## Prerequisites

Ensure you have:
- CMake 3.14+
- A C++17 compatible compiler (GCC 7.1+, Clang 5+, Intel 19+)
- FFTW and HDF5 libraries (as per ERF requirements)
- Optional: MPI for parallel testing
- Python 3.6+ (for checker scripts)

## Part 1: Build Verification

### Step 1.1: Configure the build

```bash
cd /path/to/ERF
mkdir -p build_phase11
cd build_phase11
cmake -DCMAKE_BUILD_TYPE=Release ..
```

**Expected Result**: CMake configuration succeeds with no errors about missing headers or undefined symbols.

### Step 1.2: Compile

```bash
make -j 8
```

**Expected Result**: Compilation succeeds without errors. You may see warnings about unused variables or similar (those are benign).

**Key Compilation Checks**:
- No errors in `ERF_AdvanceTwoStreamRadiation.cpp` (the main Phase 11 file)
- No undefined references to `resolve_surface_*` functions
- No missing headers (especially `<cmath>` for `std::isfinite()`)

If compilation fails, check:
1. That `#include <cmath>` is present in ERF_AdvanceTwoStreamRadiation.cpp (Line 9)
2. That all Phase 11 helper functions are correctly inlined with `AMREX_GPU_DEVICE AMREX_FORCE_INLINE`
3. That array parameter syntax in `vertical_two_stream_sweep()` matches AMReX conventions

## Part 2: Phase 11 Regression Test

### Step 2.1: Run the Phase 11 RegTest

```bash
cd /path/to/ERF/build_phase11

# Run the Phase 11 test
./erf \
  /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/inputs \
  max_level=0 \
  amr.n_cell="8 8 64" \
  stop_time=2.5 \
  erf.radiation_type="TwoStream" \
  erf.radiation.diag_file="radiation_diag_phase11.dat"
```

**Expected Runtime**: ~30 seconds to a few minutes (depends on machine and compiler)

**Expected Output**:
- Timestep progression messages (5 steps total at dt=0.5s)
- Radiation diagnostics printed to stdout at each step
- No crashes or assertion failures
- `radiation_diag_phase11.dat` file created in the run directory

### Step 2.2: Validate diagnostics output

```bash
# Check that diagnostics file was created
ls -la radiation_diag_phase11.dat

# View first 10 lines
head -10 radiation_diag_phase11.dat

# Count data rows (should be ~5, one per timestep)
grep -v "^#" radiation_diag_phase11.dat | wc -l
```

**Expected Results**:
- File size > 1000 bytes
- Header row with column names
- 5-6 data rows (one per timestep)

### Step 2.3: Run the checker script

```bash
# Make the checker executable
chmod +x /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/check_hetero_accuracy.py

# Run the checker
python3 /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/check_hetero_accuracy.py
```

**Expected Output**:
```
======================================================================
Phase 11 Surface Heterogeneity RegTest Validation
======================================================================

Working directory: /path/to/your/run/dir

1. Checking for required output files...
✓ Radiation diagnostics CSV found: ./radiation_diag_phase11.dat

2. Parsing radiation diagnostics CSV...
✓ Successfully parsed N diagnostic rows

3. Checking for NaN/Inf in diagnostics...
✓ SW_surface: all values finite
✓ SW_TOA: all values finite
✓ F_up_surface: all values finite
✓ F_down_toa: all values finite
✓ heating_rate_max: all values finite

4. Validating heating rates...
✓ Nontrivial heating detected (max = X.XXXXe-03 K/s)

5. Validating surface fluxes...
SW_surface range: [X.XX, YYY.YY] W/m²
  ✓ SW_surface in reasonable range
F_up_surface range: [X.XX, YYY.YY] W/m²
  ✓ F_up_surface in reasonable range

6. Phase 11 feature validation...
Note: Phase 11 test runs in fallback mode (hetero fields all nullptr)
- surface_albedo_sw = 0.3 (from inputs)
- surface_emissivity_lw = 0.99 (from inputs)
- surface_temp_k = 300.0 K (from inputs)
✓ Fallback path being exercised (no hetero LSM fields available)

======================================================================
VALIDATION SUMMARY: 6/6 checks passed
======================================================================

✅ ALL CHECKS PASSED - Phase 11 test successful
```

### Step 2.4: Detailed validation of diagnostics values

**Manually inspect the diagnostics CSV** for reasonableness:

```bash
cat radiation_diag_phase11.dat
```

**Expected patterns**:
- **SW_TOA**: Should equal `S0 * cos(zenith_angle)` = 1361 * cos(45°) ≈ 962 W/m²
- **SW_surface**: Should be 0 to ~600 W/m² (depends on albedo and optical depth)
  - With default albedo=0.3, expect ~(1-0.3) * incident flux
- **F_up_surface** (LW): Should be 300-500 W/m² (based on surface temp ~300K and emissivity=0.99)
  - σ*T⁴ with T=300K gives ~458 W/m²
- **heating_rate_max**: Should be 1-100 K/s (nonzero, indicating active radiation)

## Part 3: Backward Compatibility Verification

### Step 3.1: Run Phase 10 RegTest (should be unchanged)

```bash
cd /path/to/ERF/build_phase11

# Run an existing Phase 10 test (e.g., TwoStream_NonuniformDZ)
./erf \
  /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_NonuniformDZ/inputs \
  max_level=0 \
  amr.n_cell="8 8 32" \
  stop_time=2.0 \
  erf.radiation_type="TwoStream" \
  erf.radiation.diag_file="radiation_diag_phase10.dat"
```

**Expected Result**: Test completes successfully with no changes to output vs. Phase 10

### Step 3.2: Compare outputs (optional)

If you have a baseline Phase 10 run:

```bash
# Extract diagnostic values
awk '{print $4, $5}' radiation_diag_phase10.dat > phase10_fluxes.txt
awk '{print $4, $5}' radiation_diag_phase11.dat > phase11_fluxes.txt

# Compute differences (should be tiny or zero)
paste phase10_fluxes.txt phase11_fluxes.txt | \
  awk '{print $1-$3, $2-$4}' > flux_diffs.txt

head flux_diffs.txt
```

**Expected**: Differences should be < 1e-10 (floating-point rounding only)

## Part 4: GPU Safety Verification (if using GPU build)

### Step 4.1: Compile with GPU support

```bash
cd /path/to/ERF/build_phase11_gpu
cmake -DCMAKE_BUILD_TYPE=Release -DAMReX_GPU_BACKEND=CUDA ..
make -j 8
```

### Step 4.2: Run GPU test

```bash
# Run the Phase 11 test on GPU
./erf \
  /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/inputs \
  max_level=0 \
  amr.n_cell="32 32 64" \
  stop_time=1.0
```

**Expected Result**: Test completes without GPU errors (e.g., no CUDA assertion failures)

**GPU-Specific Checks**:
- No "AMREX_ALWAYS_ASSERT" failures in device code
- No "invalid device function" errors
- Results match CPU version (within floating-point tolerance)

## Part 5: Code Review Checklist

Review the implementation files to verify backward compatibility and correctness:

### ERF_RadStruct.H (check for Phase 11 fields)
- [ ] Three new fields present: `surface_albedo_sw`, `surface_emissivity_lw`, `surface_temp_k`
- [ ] Defaults are sensible: 0.3, 0.99, 300.0
- [ ] `init_params()` queries these fields from input file
- [ ] Values are clamped to valid ranges (albedo/emissivity ∈ [0,1], temp > 0)

### ERF_AdvanceTwoStreamRadiation.cpp (check Phase 11 additions)
- [ ] Five helper functions defined (all with `AMREX_GPU_DEVICE AMREX_FORCE_INLINE`):
  - `clamp_finite()`
  - `is_finite_positive()`
  - `resolve_surface_albedo_sw()`
  - `resolve_surface_emissivity_lw()`
  - `resolve_surface_temp_k()`
- [ ] `vertical_two_stream_sweep()` signature updated with 6 new parameters (3 flags, 3 arrays)
- [ ] SW surface flux uses resolved albedo (line ~730)
- [ ] LW upwelling initialization uses resolved emissivity and t_sfc (line ~640)
- [ ] `compute_twostream_radiation_diagnostics()` retrieves and passes hetero fields (line ~847-858)

### RAD_DEVELOPMENT.md (check documentation)
- [ ] Phase 11 section added before Phase 10
- [ ] Roadmap table updated (Phase 11 marked as ✅ Complete)
- [ ] Fallback chain clearly documented
- [ ] GPU safety notes present
- [ ] New RegTest name listed: `TwoStream_SurfaceHeterogeneity`

### TwoStream_SurfaceHeterogeneity RegTest (check test setup)
- [ ] Directory created: `Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/`
- [ ] Files present:
  - `inputs` (simulation configuration)
  - `input_sounding_hetero` (atmospheric sounding)
  - `check_hetero_accuracy.py` (validation script)
- [ ] Inputs file sets Phase 11 parameters and enables radiation

## Part 6: Integration Test (Optional but Recommended)

### Step 6.1: Run full radiation test suite

```bash
cd /path/to/ERF/build_phase11

# Run multiple radiation tests to verify no regressions
for test in \
  SW_ClearSky_Analytical \
  LW_Isothermal \
  SW_Cloud_Layer \
  SW_Scattering_Cloud \
  TwoStream_RhoTheta_Coupling \
  TwoStream_NonuniformDZ \
  TwoStream_SurfaceHeterogeneity; do
  
  echo "Running $test..."
  ./erf \
    /path/to/ERF/Exec/CanonicalTests/Radiation/$test/inputs \
    max_level=0 \
    stop_time=1.0 \
    2>&1 | tail -5
done
```

**Expected Result**: All tests complete successfully, no crashes or assertion failures

## Troubleshooting

### Issue: Compilation error in ERF_AdvanceTwoStreamRadiation.cpp

**Symptoms**: Undefined reference to `resolve_surface_albedo_sw` or similar

**Solution**:
1. Verify that helper functions are declared before `vertical_two_stream_sweep()`
2. Check that `AMREX_GPU_DEVICE AMREX_FORCE_INLINE` is present for all helpers
3. Ensure `std::isfinite()` is available (check for `#include <cmath>`)

### Issue: Runtime error: "Array access out of bounds"

**Symptoms**: Segfault or assertion failure during device kernel

**Solution**:
1. Verify that `Array4::contains()` checks are used before accessing hetero arrays
2. Check that nullptr checks are done on `hetero_alb_sw`, `hetero_emiss_lw`, `t_sfc`
3. Ensure device function parameters are correct (const arrays, proper lambda captures)

### Issue: Diagnostics values are NaN or Inf

**Symptoms**: CSV file contains NaN or very large values

**Solution**:
1. Check that `clamp_finite()` is being called on hetero field values
2. Verify that fallback values are always finite (defaults are hardcoded)
3. Ensure surface temperature is positive (checked by `is_finite_positive()`)

## Success Criteria

✅ All checks pass when:

1. **Compilation**: No errors or warnings in Phase 11 code
2. **Runtime**: Phase 11 RegTest completes without crashes
3. **Diagnostics**: CSV file contains valid, finite flux and heating values
4. **Validation**: Checker script reports "6/6 checks passed"
5. **Backward Compat**: Phase 10 and earlier tests produce identical outputs
6. **Physical Realism**: Heating rates are nonzero, surface fluxes in expected ranges

## Expected Diagnostic Values

For the default Phase 11 test configuration:

```
Config:
- Domain: 3000m × 3000m × 1024m
- Grid: 8 × 8 × 64 cells
- Zenith angle: 45°
- S0 = 1361 W/m²
- tau_sw = 0.05 per layer
- tau_lw = 1.0 per layer
- surface_albedo_sw = 0.3
- surface_emissivity_lw = 0.99
- surface_temp_k = 300.0 K

Expected Outputs:
- SW_TOA ≈ 962 W/m² (1361 * cos(45°))
- SW_surface ≈ 300-400 W/m² (depends on optical depth)
- F_up_surface ≈ 450 W/m² (≈0.99 * σ * 300^4)
- heating_rate_max ≈ 1-10 K/s (typical for clear-sky)
```

## Next Steps (Future Phases)

Phase 11 implementation enables future work:

- **Phase 12**: Dynamically diagnosed τ(k) from moisture and cloud fields
- **Phase 13**: PBL coupling with spatially varying albedo/emissivity
- **Phase 14**: Prognostic surface temperature evolution
- **Phase 15+**: Aerosol effects, advanced microphysics coupling

## Questions or Issues?

If you encounter problems:

1. Check the compiler output for specific error messages
2. Verify that all modified files have been committed
3. Compare your local files with the git repository
4. Review the Phase 11 section in `RAD_DEVELOPMENT.md` for implementation details
5. Examine the checker script output for clues about what failed

---

**Document Version**: 1.0 (Phase 11 Initial Release)  
**Date**: 2026-08-08  
**Maintainer**: ERF Radiation Development Team
