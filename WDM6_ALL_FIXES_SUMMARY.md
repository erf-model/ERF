# WDM6 Complete Fix Summary

**Date**: 2026-08-02  
**Status**: Ready for testing  
**Branch**: HEAD

---

## Changes Applied

### 1. **CRITICAL FIX: nn Initialization** ✅ FIXED

**File**: `Source/Microphysics/WDM6/ERF_InitWDM6.cpp` (Line 69)

**Problem**: nn was density-dependent, causing 12x higher values at altitude  
**Impact**: Runaway nc accumulation → tiny droplets → storm collapse by hour 2

**Change:**
```cpp
// Before (WRONG):
nn(i,j,k) = ccn0_init / rho;  // Varies 15x from surface to upper troposphere!

// After (CORRECT):
nn(i,j,k) = ccn0_init;  // Constant everywhere, matches WRF
```

**Expected Result**: Rain persists through hour 2, cold pool maintained

---

### 2. **Ice Species Initialization** ✅ RESTORED

**File**: `Source/Initialization/ERF_InitFromWRFInput.cpp` (Lines 247-262)

**Change**: Restored reading QICE, QSNOW, QGRAUP from wrfinput for WDM6

**Before:**
```cpp
if (n_qstate_moist >= 11) {
    // Only Morrison reads ice species
    NC_names.push_back("QICE");
    NC_names.push_back("QRAIN");
    NC_names.push_back("QSNOW");
    NC_names.push_back("QGRAUP");
} else if (n_qstate_moist >= 6) {
    // WDM6/WSM6: skip ice species, only read QRAIN
    NC_names.push_back("QRAIN");
}
```

**After:**
```cpp
if (n_qstate_moist >= 6) {
    // All 6+ class schemes (WDM6, WSM6, Morrison): read all ice species
    NC_names.push_back("QICE");
    NC_names.push_back("QRAIN");
    NC_names.push_back("QSNOW");
    NC_names.push_back("QGRAUP");
}
```

**Why**: Original optimization (skip ice init, let physics generate) caused memory issues. Restoring full initialization from wrfinput.

---

### 3. **Previous Fixes** (Already in codebase)

These were applied in earlier commits and are already working:

#### 3a. **Theta Update Fix** (Commit 5a692b12)
- **File**: `Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp` (Lines 439-453)
- **Issue**: Theta wasn't being updated after Fortran microphysics
- **Fix**: Convert T → theta after Fortran call, update state arrays
- **Status**: ✅ Already fixed

#### 3b. **Storm Cell Diagnostics** (Commit 1df0727c)
- **File**: `Source/Microphysics/WDM6/ERF_UpdateWDM6.cpp`
- **Added**: Print diagnostics at storm core (i=240, j=150)
- **Shows**: T, P, den, qv, qc, qr, qi, nn, nc, nr before/after
- **Status**: ✅ Already working (used for debugging!)

#### 3c. **nn State Management** (Commit 1df0727c)
- **File**: `Source/Microphysics/WDM6/ERF_InitWDM6.cpp`
- **Issue**: nn was being overwritten by zero from state on first call
- **Fix**: Initialize nn in Init(), skip state read on first Copy_State_to_Micro
- **Status**: ✅ Already fixed (but initialization value was wrong - fixed today!)

---

## Files Modified Today (2026-08-02)

1. **Source/Microphysics/WDM6/ERF_InitWDM6.cpp**
   - Line 69: Removed `/rho` from nn initialization
   - Added comment explaining WRF convention

2. **Source/Initialization/ERF_InitFromWRFInput.cpp**
   - Lines 247-262: Restored ice species reading for WDM6
   - Lines 526-562: Already had correct QICE/QSNOW/QGRAUP mapping

---

## Build and Test Instructions

### Step 1: Rebuild ERF
```bash
cd /g/g10/wise14/compiling/clean/ERF
cmake --build build -j8
```

### Step 2: Run Test Case
Use the same deep convective storm case that showed the problem.

### Step 3: Check Diagnostics

**At initialization (step 1):**
```bash
grep "nn (#/kg)" <erf_log> | head -1
```
**Expected**: `nn max ≈ 1.0e8` (NOT 1.2e9!)

**At hour 1 (step 3600):**
```bash
grep "qr (kg/kg)" <erf_log> | grep "Timestep 3600"
```
**Expected**: `qr max ≈ 7-8 g/kg` (NOT 4 g/kg!)

**At hour 2 (step 7200):**
```bash
grep "qr (kg/kg)" <erf_log> | grep "Timestep 7200"
```
**Expected**: `qr max > 0.5 g/kg` (NOT 0!)

**Storm cell diagnostics:**
```bash
grep "Storm cell diagnostics" <erf_log> | head -5
```
Should show reasonable nn/nc/nr values matching WRF.

### Step 4: Visual Verification

Plot theta at surface at hour 2:
- **Should see**: Cold pool (5-10K cooler where rain is falling)
- **Should NOT see**: Uniform warm temperature (no cold pool)

---

## Expected Behavior After Fixes

### Comparison Table (Hour 1 - Step 3600):

| Metric | WRF | ERF Before | ERF After (Expected) |
|--------|-----|------------|----------------------|
| **nn init** | 1.0e8 | **1.2e9** ❌ | **1.0e8** ✓ |
| **qr max** | 7.9 g/kg | **4.1 g/kg** ⚠️ | **~7.5 g/kg** ✓ |
| **nc max** | 4.2e10 | 3.4e10 | **~4.0e10** ✓ |

### Comparison Table (Hour 2 - Step 7200):

| Metric | WRF | ERF Before | ERF After (Expected) |
|--------|-----|------------|----------------------|
| **qr max** | 0.70 g/kg ✓ | **0** ❌ | **~0.6 g/kg** ✓ |
| **Cold pool** | Present ✓ | **GONE** ❌ | **Present** ✓ |
| **Surface T** | Cool ✓ | **Too warm** ❌ | **Cool** ✓ |

---

## What Was the Root Cause?

### The Bug Chain:

1. **ERF initialized nn wrong**: `nn = CCN0/rho` → 12x too high at altitude
2. **High nn at altitude** → Very high nc activation in upper-level clouds
3. **Very high nc** → Droplets too small to grow efficiently
4. **Small droplets** → Weak autoconversion, less rain production
5. **Less rain at hour 1** → Only 50% of WRF's rain mass
6. **Small rain amount** → Completely evaporates/sediments by hour 2
7. **No rain** → No evaporative cooling → Cold pool collapses
8. **No cold pool** → Storm dies

### The Fix:

Simply match WRF's initialization: `nn = CCN0` (constant with height)

This single-line change breaks the runaway feedback loop!

---

## Additional Documentation

Created today:
1. **WDM6_RUNAWAY_NN_BUG_FIX.md** - Detailed analysis with diagnostics
2. **WDM6_CPU_FORTRAN_DIAGNOSTICS.md** - Diagnostic procedures
3. **WDM6_ALL_FIXES_SUMMARY.md** - This file

Previous documentation (from earlier work):
1. **WDM6_COMPREHENSIVE_DIAGNOSTICS.md**
2. **WDM6_COUPLING_VERIFICATION.md**
3. **WDM6_DEFAULT_ISSUES.md**
4. **WDM6_DIAGNOSTICS_FINAL.md**
5. **WDM6_DIAGNOSTIC_COMPARISON.md**
6. **WDM6_FIXES_APPLIED.md**
7. **WDM6_VS_MORRISON_COMPARISON.md**

---

## Commit Message (Suggested)

```
[WDM6] Fix critical nn initialization causing storm collapse

Root cause: nn was initialized as ccn0/rho, varying 15x with altitude.
At upper levels (rho=0.083), nn was 12x higher than WRF, causing runaway
nc accumulation, tiny droplets, and complete precipitation loss by hour 2.

Fix: Initialize nn as constant (matching WRF), breaking feedback loop.

Also restored QICE/QSNOW/QGRAUP initialization from wrfinput for WDM6.

Files changed:
- Source/Microphysics/WDM6/ERF_InitWDM6.cpp (line 69)
- Source/Initialization/ERF_InitFromWRFInput.cpp (lines 247-262)

Expected result: Rain and cold pool persist through hour 2, matching WRF.

Co-Authored-By: Claude <noreply@anthropic.com>
```

---

## Next Steps

1. ✅ **Code changes applied** (today)
2. ⏳ **Rebuild ERF** (user action needed)
3. ⏳ **Rerun test case** (user action needed)
4. ⏳ **Verify diagnostics** (user action needed)
5. ⏳ **Commit changes** (if tests pass)
6. ⏳ **Run longer simulations** (6+ hours) to confirm stability

---

**Status**: Code ready, awaiting user testing

**Confidence**: 95% - The diagnostic evidence strongly supports this fix

**Risk**: Low - Changes are minimal and match WRF behavior exactly

---

**Generated**: 2026-08-02  
**By**: Claude Code (Sonnet 4)
