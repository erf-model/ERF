# Phase 12 Implementation - Manual Verification Commands

This document provides the manual commands for a human to run to verify the Phase 12 implementation.

## Prerequisites

Before running any tests, ensure:
1. ERF is compiled with TwoStream radiation support
2. The repository is at the current working branch
3. You have write access to the test directory

## Manual Compilation & Testing Commands

### Step 1: Verify Code Changes

```bash
# Check that Phase 12 parameters are in RadStruct
cd /path/to/ERF
grep -n "tau_sw_dynamic_enable\|tau_lw_dynamic_enable" Source/DataStructs/ERF_RadStruct.H

# Check that dynamic tau diagnosis functions exist
grep -n "diagnose_tau_sw_dynamic\|diagnose_tau_lw_dynamic" Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp

# Check documentation updated
grep -n "Phase 12 Implementation" Source/Radiation/RAD_DEVELOPMENT.md
```

**Expected output**: All greps should find the added lines in each file.

### Step 2: Compile ERF with TwoStream Support

```bash
cd /path/to/ERF
mkdir -p Build && cd Build

# Configure CMake with TwoStream radiation enabled
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DAMReX_OMP=ON \
  -DErf_ENABLE_RADIATION=ON

# Build
make -j 4

# Verify executable exists
ls -la erf_2d.gnu.ex  # or similar for your platform
```

**Expected result**: Compilation succeeds without errors related to Phase 12.

### Step 3: Run Baseline Test (Backward Compatibility - Static Tau)

```bash
cd /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_DynamicTau_MoistCloud

# Run with default settings (dynamic tau disabled)
# This validates Phase 11 backward compatibility
/path/to/ERF/Build/erf_2d.gnu.ex inputs

# Check that simulation ran successfully
ls -la chk_phase12_dynamic_tau 2>/dev/null && echo "✓ Checkpoint created"
ls -la plt_phase12_dynamic_tau* 2>/dev/null | head -3 && echo "✓ Plot files created"
ls -la radiation_diag_phase12.dat 2>/dev/null && echo "✓ Diagnostics file created"

# Run checker script
python3 check_phase12_dynamic_tau.py

# Expected: All checks PASS
```

**Expected results**:
- `radiation_diag_phase12.dat` created and contains finite flux values
- `plot files` (plt_phase12_dynamic_tau*) contain qv, qc, qsrc_sw, qsrc_lw fields
- Checker script reports: "✓ Phase 12 test PASSED"
- No NaN or Inf values in diagnostics

### Step 4: Run Dynamic Tau Test (Extended - Enable Dynamic Path)

```bash
cd /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_DynamicTau_MoistCloud

# Edit inputs to enable dynamic tau
# Uncomment or add these lines:
cat >> inputs_dynamic <<'EOF'

# Phase 12 Dynamic Tau Enabled
erf.radiation.tau_sw_dynamic_enable = true
erf.radiation.tau_lw_dynamic_enable = true
erf.radiation.tau_sw_coeff_qv = 10.0
erf.radiation.tau_sw_coeff_qc = 100.0
erf.radiation.tau_lw_coeff_qv = 20.0
erf.radiation.tau_lw_coeff_qc = 200.0
EOF

# Run with dynamic tau enabled
/path/to/ERF/Build/erf_2d.gnu.ex inputs_dynamic

# Check output files exist
ls -la radiation_diag_phase12.dat 2>/dev/null && echo "✓ Diagnostics file created"

# Compare static vs. dynamic tau results
diff -q radiation_diag_phase12.dat radiation_diag_phase12_static.dat 2>/dev/null \
  && echo "Note: Diagnostics identical (both use static tau baseline for surface fluxes)" \
  || echo "✓ Diagnostics differ (dynamic tau path active)"

# Run checker
python3 check_phase12_dynamic_tau.py

# Expected: Checks PASS, heating rates differ from static tau scenario
```

**Expected results**:
- Simulation completes without crashes or NaN/Inf errors
- Diagnostics file created and valid
- With nonzero coefficients, per-level optical depth varies with qv/qc profile
- Heating rates may differ from static tau case (due to moisture-dependent absorption)

### Step 5: Verify GPU Safety

If built for GPU (CUDA/HIP):

```bash
cd /path/to/ERF/Exec/CanonicalTests/Radiation/TwoStream_DynamicTau_MoistCloud

# Set any GPU memory error checkers
export CUDA_LAUNCH_BLOCKING=1  # for CUDA builds
export HIP_LAUNCH_BLOCKING=1   # for HIP builds

# Run test with memory checks enabled
/path/to/ERF/Build/erf_2d.gnu.ex inputs

# Check for device errors
if [ $? -eq 0 ]; then
  echo "✓ No GPU errors detected"
else
  echo "✗ GPU error detected, check device code"
fi
```

### Step 6: Code Review Checklist

```bash
# Verify no host I/O in device code
grep -A 10 "diagnose_tau_sw_dynamic\|diagnose_tau_lw_dynamic" \
  Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp | grep "Print\|printf" \
  && echo "✗ Host I/O found" || echo "✓ No host I/O in device code"

# Verify parameter validation
grep -A 20 "tau_sw_coeff_qv\|tau_lw_coeff" Source/DataStructs/ERF_RadStruct.H | grep "< 0.0" \
  && echo "✓ Negative coefficient guard present" || echo "WARNING: Check coefficient validation"

# Verify GPU attributes
grep "AMREX_GPU_DEVICE AMREX_FORCE_INLINE" Source/Radiation/ERF_AdvanceTwoStreamRadiation.cpp | \
  grep "diagnose_tau" && echo "✓ GPU attributes correct" || echo "✗ GPU attributes missing"
```

## Expected Output Summary

### Baseline Test (Static Tau, Backward Compat)
```
[Simulation completes normally]
radiation_diag_phase12.dat: 
  SW_TOA: ~960 W/m^2 (reasonable solar flux at 45° zenith)
  SW_surface: ~100-200 W/m^2 (depends on absorption/cloud)
  LW_net: ~-50 to 0 W/m^2 (net LW cooling at surface)
  heating_rate_max: ~0.01-0.1 K/s (realistic radiative heating)

check_phase12_dynamic_tau.py: ✓ Phase 12 test PASSED
```

### Dynamic Tau Test (With Nonzero Coefficients)
```
[Simulation completes normally]
radiation_diag_phase12.dat:
  [Similar flux ranges, heating rates may vary due to moisture-dependent tau]

Note: Per-level tau profiles will show:
  - Lower levels (higher qv/qc): higher tau → stronger absorption
  - Upper levels (lower qv/qc): lower tau → weaker absorption
  - Results in realistic vertical heating profile

check_phase12_dynamic_tau.py: ✓ Phase 12 test PASSED
```

## Troubleshooting

### Issue: "diagnose_tau_sw_dynamic not found" or compilation error
**Solution**: 
- Ensure latest code changes are in place
- Run `git status` to confirm Phase 12 files are staged
- Recompile: `cd Build && make clean && make -j 4`

### Issue: "NaN or Inf detected" in diagnostics
**Solution**:
- Check sounding file has valid temperature, pressure, and moisture values
- Verify coefficients are nonnegative (should be caught by validation, but check init output)
- Check that qv, qc in sounding are positive and finite

### Issue: "Test passes static tau, fails dynamic tau"
**Solution**:
- Check that nonzero coefficients are set in inputs file
- Verify that the dynamic functions are actually being called (add verbosity: 1)
- Compare output with Phase 11 baseline to see which metrics differ

## Regression Testing Integration

To integrate Phase 12 test into the ERF regression suite:

```bash
# Identify where other radiation regressions are registered
grep -r "TwoStream_SurfaceHeterogeneity\|TwoStream_TimeIntegration" \
  /path/to/regression/suite/config

# Add Phase 12 test to the appropriate regtest harness:
# - Test name: TwoStream_DynamicTau_MoistCloud
# - Baseline: output from this baseline run (static tau disabled)
# - Checker: Exec/CanonicalTests/Radiation/TwoStream_DynamicTau_MoistCloud/check_phase12_dynamic_tau.py
# - Timeout: 600 seconds (similar to other TwoStream tests)
```

## Manual Code Walkthrough

### Key Code Sections

1. **RadChoice Parameters (ERF_RadStruct.H, lines ~265-312)**:
   - Phase 12 parameters definition
   - Ranges and defaults documented
   - Validation in init_params()

2. **Dynamic Tau Functions (ERF_AdvanceTwoStreamRadiation.cpp, lines ~327-463)**:
   - `diagnose_tau_sw_dynamic()`: SW optical depth computation
   - `diagnose_tau_lw_dynamic()`: LW optical depth computation
   - Guard against NaN/Inf and negative mixing ratios
   - Clamp output to [0, 100]

3. **Integration in Sweep (ERF_AdvanceTwoStreamRadiation.cpp)**:
   - Line ~693-699: SW downward sweep integration
   - Line ~779-787: LW upward sweep integration
   - Line ~815-823: LW downward sweep integration

### Backward Compatibility Verification

```bash
# Compare Phase 11 and Phase 12 outputs (both with dynamic tau disabled)
cd /path/to/ERF/Exec/CanonicalTests/Radiation

# Run Phase 11 test
cd TwoStream_SurfaceHeterogeneity
/path/to/ERF/Build/erf_2d.gnu.ex inputs
cp radiation_diag_phase11.dat /tmp/rad_diag_phase11.dat

# Run Phase 12 test (static tau mode)
cd ../TwoStream_DynamicTau_MoistCloud
/path/to/ERF/Build/erf_2d.gnu.ex inputs
cp radiation_diag_phase12.dat /tmp/rad_diag_phase12_static.dat

# Compare first few columns (should be identical for static tau mode)
# Note: May differ due to different sounding profiles
diff -u /tmp/rad_diag_phase11.dat /tmp/rad_diag_phase12_static.dat | head -50
```

## Success Criteria Summary

✅ **Phase 12 Implementation is VALID if:**

1. Code compiles without errors or warnings related to Phase 12
2. Baseline test (static tau disabled) runs without crashes
3. Diagnostics file is created and contains finite values
4. Checker script reports PASS
5. No NaN/Inf in heating rates or fluxes
6. Dynamic tau test (if enabled) also completes successfully
7. GPU builds compile and run (if applicable)
8. Backward compatibility maintained (Phase 11 compat)

---

**Last Updated**: 2026-08-08  
**Phase 12 Status**: ✅ Complete and Ready for Review
