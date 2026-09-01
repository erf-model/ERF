# Phase 11 Implementation - Human Verification Instructions

This document provides **exact commands** for a human to verify the Phase 11 implementation on their local system. **The agent has NOT executed these commands** — this is for you to run.

## Quick Start Verification (5 minutes)

### 1. Verify Code Compilation

```bash
cd /path/to/ERF
rm -rf build_phase11_verify
mkdir build_phase11_verify && cd build_phase11_verify
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 8 2>&1 | grep -E "error:|undefined|missing" || echo "✓ Compilation successful"
```

**Expected**: No error messages, build completes.

### 2. Verify Phase 11 Test Runs

```bash
cd /path/to/ERF/build_phase11_verify

./erf \
  ../Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/inputs \
  max_level=0 \
  amr.n_cell="8 8 64" \
  stop_time=2.5 \
  2>&1 | tail -20
```

**Expected**: Test completes, no crashes. Look for "AMR: Grid created" and timestep progression messages.

### 3. Verify Diagnostics File Created

```bash
cd /path/to/your/run/directory
ls -la radiation_diag_phase11.dat
head -20 radiation_diag_phase11.dat
```

**Expected**: File exists, contains header and data rows.

### 4. Run Validation Script

```bash
cd /path/to/your/run/directory
python3 /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/check_hetero_accuracy.py
```

**Expected Output**:
```
VALIDATION SUMMARY: 6/6 checks passed
✅ ALL CHECKS PASSED - Phase 11 test successful
```

## Full Verification (30-60 minutes)

### Step 1: Build and Run Phase 11 Test

```bash
cd /path/to/ERF
rm -rf build_test
mkdir build_test && cd build_test
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 8

# Run Phase 11 RegTest
./erf \
  ../Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/inputs \
  max_level=0 \
  stop_time=2.5 \
  2>&1 | tee phase11_run.log

# Check for success
tail -50 phase11_run.log | grep -E "Success|PASSED" || echo "Check run.log for completion"
```

**Expected**: Run completes in 30 seconds to 2 minutes.

### Step 2: Validate Phase 11 Output

```bash
# Check diagnostics file
ls -la radiation_diag_phase11.dat
wc -l radiation_diag_phase11.dat  # Should have 6-8 lines (header + 5-6 data rows)

# Inspect CSV data
head -20 radiation_diag_phase11.dat

# Run checker script
cd /path/to/run/directory
python3 /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/check_hetero_accuracy.py

# Verify exact values
awk -F, 'NR>2 {print $4, $5, $6, $7, $8}' radiation_diag_phase11.dat
```

**Expected SW_TOA (~962), SW_surface (~300-400), F_up (~450), heating_max (~1-10)**

### Step 3: Test Backward Compatibility

```bash
# Run Phase 10 test (should be identical to Phase 11 output when hetero=nullptr)
cd /path/to/ERF/build_test

./erf \
  ../Exec/CanonicalTests/Radiation/TwoStream_NonuniformDZ/inputs \
  max_level=0 \
  stop_time=2.0 \
  2>&1 | tee phase10_run.log

# Compare diagnostics
tail -10 phase10_run.log | grep "heating_rate_max\|SW_surface"
```

**Expected**: Output should match Phase 10 baseline (within floating-point tolerance).

### Step 4: Code Review Spot Checks

```bash
# Verify Phase 11 parameters in RadStruct
grep -A 5 "surface_albedo_sw" /path/to/ERF/Source/DataStructs/ERF_RadStruct.H

# Verify helper functions exist
grep "resolve_surface_" /path/to/ERF/Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp | head -5

# Verify GPU annotations
grep -c "AMREX_GPU_DEVICE" /path/to/ERF/Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp
# Expected: at least 10 (5+ new helpers + existing functions)

# Verify documentation
grep -A 10 "Phase 11" /path/to/ERF/Source/Radiation/RAD_DEVELOPMENT.md | head -15
```

**Expected**: All references found, GPU annotations present.

## Extended Verification (60+ minutes)

### Run All Radiation Tests

```bash
cd /path/to/ERF/build_test

tests=(
  "SW_ClearSky_Analytical"
  "LW_Isothermal"
  "SW_Cloud_Layer"
  "SW_Scattering_Cloud"
  "TwoStream_RhoTheta_Coupling"
  "TwoStream_NonuniformDZ"
  "TwoStream_SurfaceHeterogeneity"
)

for test in "${tests[@]}"; do
  echo "===== Testing $test ====="
  ./erf \
    ../Exec/CanonicalTests/Radiation/$test/inputs \
    max_level=0 \
    stop_time=1.0 \
    amr.n_cell="8 8 32" \
    2>&1 | tail -5 | grep -E "Success|PASSED|Error" || echo "Test $test completed"
  echo ""
done
```

**Expected**: All tests complete successfully, no errors.

### GPU Build Verification (if CUDA/GPU available)

```bash
cd /path/to/ERF
rm -rf build_gpu
mkdir build_gpu && cd build_gpu

cmake -DCMAKE_BUILD_TYPE=Release -DAMReX_GPU_BACKEND=CUDA ..
make -j 8 2>&1 | grep -i "cuda\|gpu"

# Run Phase 11 test on GPU
./erf \
  ../Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/inputs \
  max_level=0 \
  amr.n_cell="16 16 64" \
  stop_time=1.0 \
  2>&1 | tee gpu_run.log

# Check for GPU errors
grep -i "error\|assertion" gpu_run.log || echo "✓ No GPU errors detected"
```

**Expected**: GPU build succeeds, test runs without GPU errors.

## Diagnostic Commands for Troubleshooting

### If test crashes:

```bash
# Get full error traceback
./erf ../Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/inputs 2>&1 | head -100

# Check for specific Phase 11 error messages
grep -i "phase.11\|surface\|albedo\|emiss" run.log
```

### If diagnostics file missing or empty:

```bash
# Check current directory
pwd
ls -la radiation_diag*.dat

# Look for the file in various locations
find . -name "radiation_diag*.dat" 2>/dev/null

# Check file permissions
stat radiation_diag_phase11.dat 2>/dev/null
```

### If checker script fails:

```bash
# Debug CSV parsing
head -50 radiation_diag_phase11.dat | tail -20

# Check for NaN/Inf manually
awk -F, 'NR>2 {for(i=4; i<=NF; i++) if ($i ~ /NaN|Inf/) print "Row " NR ": " $i}' radiation_diag_phase11.dat

# Get summary stats
awk -F, 'NR>2 {print "SW_surface: " $4 ", heating: " $8}' radiation_diag_phase11.dat
```

## Documentation to Review

After successful verification, review these documents:

1. **Source/Radiation/RAD_DEVELOPMENT.md** - Read Phase 11 section (~150 lines)
   - Understand fallback chain
   - Review GPU safety measures
   - Check backward compatibility notes

2. **Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/README.md** - Test overview
   - Test configuration and purpose
   - Expected output values
   - Future extension notes

3. **Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/PHASE_11_SUMMARY.md** - Complete summary
   - All files modified (with line counts)
   - Code quality checks
   - Implementation checklist

4. **Exec/CanonicalTests/Radiation/TwoStream_SurfaceHeterogeneity/MANUAL_VERIFICATION.md** - Detailed procedures
   - Step-by-step verification
   - Troubleshooting guide
   - Success criteria

## Final Checklist

Use this checklist to confirm Phase 11 is correctly implemented:

- [ ] Code compiles without errors
- [ ] Phase 11 RegTest runs to completion
- [ ] Diagnostics file is created with valid data
- [ ] All checker script checks pass (6/6)
- [ ] SW_TOA ≈ 962 W/m² (±5%)
- [ ] SW_surface ≈ 300-400 W/m² (±10%)
- [ ] F_up_surface ≈ 450 W/m² (±5%)
- [ ] heating_rate_max ≈ 1-10 K/s (nonzero)
- [ ] Phase 10 tests show identical output
- [ ] All documentation files present and readable
- [ ] GPU tests pass (if GPU build available)

## Expected Timings

- Quick verification: ~5 minutes
- Full verification: ~30-60 minutes
- Extended verification (all tests + GPU): ~1-2 hours

## Success Criteria

Phase 11 is **successfully verified** when:

✅ All commands execute without errors  
✅ Phase 11 RegTest completes successfully  
✅ Diagnostics contain physically sensible values  
✅ Checker script reports "6/6 checks passed"  
✅ All documentation is present and coherent  
✅ Previous tests (Phase 1-10) show no regressions  

## Next Steps

Once verification is complete:

1. Create a Pull Request with the Phase 11 changes
2. Note that this is Phase 11: Surface Heterogeneity + Fallback
3. Include link to this verification guide
4. Refer to RAD_DEVELOPMENT.md for technical details
5. Plan Phase 12: Dynamic Optical Depth (next phase)

---

**Created**: 2026-08-08  
**For**: Phase 11 Implementation Verification  
**Status**: Ready for Human Execution
