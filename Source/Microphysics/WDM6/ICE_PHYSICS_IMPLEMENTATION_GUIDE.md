# WDM6 Ice Physics Implementation Guide

**Status:** Warm rain GPU physics complete. Ice physics needs to be ported from WSM6.

**Goal:** Add ice processes (qi, qs, qg) to complete GPU WDM6 implementation.

**Estimated time:** 2-3 hours

**Difficulty:** Medium (copy/paste from WSM6 with minimal changes)

---

## Current State

### ✅ What's Done
- Warm rain physics (qc, qr with nc, nr number concentrations)
- CCN activation (nn → nc)
- Double-moment autoconversion
- Number-coupled rain sedimentation
- All helper functions for warm rain
- File compiles and runs on GPU

### ❌ What's Missing
- Ice depositional growth (Bergeron, vapor deposition)
- Freezing processes (homogeneous, heterogeneous, immersion)
- Melting processes (snow→rain, ice→rain, graupel→rain)
- Ice-ice interactions (aggregation, riming)
- Ice sedimentation (qi, qs, qg terminal velocities)

---

## Implementation Strategy

**Key Insight:** Ice physics is IDENTICAL in WSM6 and WDM6! Both use single-moment for ice.
Only warm rain is double-moment in WDM6 (which we already implemented).

Therefore: **Copy WSM6 ice code verbatim** (~1400 lines)

---

## Step-by-Step Instructions

### Step 0: Preparation

```bash
cd /g/g10/wise14/compiling/clean/ERF
git checkout WDM6
git pull
cd Source/Microphysics/WDM6
```

**Files you'll edit:**
- `ERF_AdvanceWDM6.cpp` (currently 650 lines → will become ~2050 lines)

**Reference file:**
- `../WSM6/ERF_AdvanceWSM6.cpp` (2508 lines)

---

### Step 1: Add Saturation Helper Functions (15 min)

**Location in target:** After line 70 (after existing helper functions)

**Copy from WSM6:** Lines 72-165

**What to copy:**
```cpp
// Saturation vapor pressure functions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real qvs_water(...) { ... }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real qvs_ice(...) { ... }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real fpvs(...) { ... }
```

**Action:**
1. Open `../WSM6/ERF_AdvanceWSM6.cpp`
2. Copy lines 72-165 (saturation functions)
3. Paste after line 70 in `ERF_AdvanceWDM6.cpp`
4. No modifications needed - use as-is

**Search for in WSM6:** `"AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE"` around line 72
**Look for:** `qvs_water`, `qvs_ice`, `fpvs` functions

---

### Step 2: Add Ice Terminal Velocity Functions (10 min)

**Location in target:** After saturation functions

**Copy from WSM6:** Lines 166-250

**What to copy:**
```cpp
// Terminal velocity calculations
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real fall_qi(...) { ... }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real fall_qs(...) { ... }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real fall_qg(...) { ... }
```

**Action:**
1. Copy WSM6 lines 166-250
2. Paste after saturation functions
3. No modifications needed

**Search for in WSM6:** `"fall_qi"`, `"fall_qs"`, `"fall_qg"`

---

### Step 3: Add Ice Process Kernels (Main Work - 60-90 min)

**Location in target:** After line 595 (after "Step 9: Evaporation", before "Step 10: Sedimentation")

**Copy from WSM6:** Lines 1400-2400 (~1000 lines)

This is the bulk of the work. Copy these sections in order:

#### 3A: Depositional Growth (Lines 1400-1600)

**Search in WSM6 for:** `"Step X: Depositional growth"` or `"Bergeron"`

**What it does:**
- Ice grows at expense of cloud water (Bergeron process)
- Vapor deposition on snow
- Vapor deposition on graupel

**Copy the entire ParallelFor block:**
```cpp
// Depositional growth
ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    // ... ice growth code ...
});
```

**Insert location:** After line 595 in WDM6

#### 3B: Freezing Processes (Lines 1600-1800)

**Search in WSM6 for:** `"freezing"` or `"piacr"` (rain freezing)

**What it does:**
- Homogeneous freezing (T < -40°C)
- Heterogeneous/immersion freezing
- Rain → graupel freezing

**Copy the ParallelFor block**

#### 3C: Melting Processes (Lines 1800-2000)

**Search in WSM6 for:** `"melting"` or `"psmlt"` (snow melt)

**What it does:**
- Snow → rain (below 0°C)
- Ice → rain
- Graupel → rain

**Copy the ParallelFor block**

#### 3D: Ice-Ice Interactions (Lines 2000-2200)

**Search in WSM6 for:** `"aggregation"` or `"riming"` or `"pracs"` (rain-snow)

**What it does:**
- Snow aggregation (snow + cloud → snow)
- Riming (cloud + graupel)
- Rain collection by snow/graupel

**Copy the ParallelFor block**

#### 3E: Ice Sedimentation (Lines 2200-2400)

**Search in WSM6 for:** `"ice sedimentation"` or `"fall_qi"`

**What it does:**
- Terminal velocity calculations
- Top-down fall for qi, qs, qg
- Snow and graupel accumulation

**Copy the ParallelFor block**

**Important:** This replaces the simplified rain sedimentation currently at line 597-637.
- Delete lines 597-637 (simplified sedimentation placeholder)
- Replace with full sedimentation from WSM6 that handles qi, qr, qs, qg

---

### Step 4: Update Step Numbers (5 min)

After inserting ice processes, you'll have many "Step X" comments that are out of order.

**Action:** Renumber the step comments sequentially:
- Step 10 becomes the first ice process
- Continue numbering through all ice processes
- Final sedimentation becomes the last step number

**Optional but recommended** - helps with debugging later.

---

### Step 5: Verify Array Access (10 min)

**Check that these arrays are used correctly:**

Already defined in WDM6:
- ✅ `qi_arr` (ice)
- ✅ `qs_arr` (snow)  
- ✅ `qg_arr` (graupel)
- ✅ `snow_arr` (snow accumulation)
- ✅ `graup_arr` (graupel accumulation)

**Action:** Search copied code for these array names and verify they match the WDM6 names.

Most should be identical to WSM6, but double-check:
```bash
grep -n "qi_arr\|qs_arr\|qg_arr" ERF_AdvanceWDM6.cpp | head -20
```

Should see them used in the ice kernels you just copied.

---

### Step 6: Build and Test (15 min)

```bash
cd /g/g10/wise14/compiling/clean/ERF

# Clean build
rm -rf build

# Configure for GPU (no Fortran)
cmake -DERF_ENABLE_WDM6_FORT=OFF \
      -DERF_PRECISION=DOUBLE \
      -B build

# Build
cmake --build build -j16 2>&1 | tee build.log

# Check for errors
tail -100 build.log
```

**Expected result:** Should compile successfully

**Common errors and fixes:**

**Error:** `undefined reference to 'some_function'`
- **Fix:** Missing helper function - copy from WSM6

**Error:** `'variable' was not declared in this scope`
- **Fix:** Variable name mismatch - check WSM6 vs WDM6 naming

**Error:** `no matching function for call to 'ParallelFor'`
- **Fix:** Box dimension mismatch - use `box` for 3D, `amrex::makeSlab(box,2,klo)` for 2D

---

### Step 7: Validate Results (Optional - 30 min)

**Run a test case with ice:**

```bash
cd Tests/WDM6_verification
./build/erf3d inputs_supercell  # or any ice-containing test
```

**What to check:**
1. qi, qs, qg values are > 0 in cloud regions
2. No NaN values
3. Snow and graupel accumulation at surface
4. Compare with Fortran version (if available)

**Basic validation:**
```bash
# Check plotfile for reasonable values
python analyze_wdm6.py plt_supercell00100

# Should show:
# qi: 1e-6 to 1e-3 kg/kg (ice mixing ratio)
# qs: 1e-5 to 1e-2 kg/kg (snow mixing ratio)  
# qg: 1e-5 to 1e-2 kg/kg (graupel mixing ratio)
# snow_accum: 0 to 100 mm (accumulated snow)
```

---

## Detailed Line-by-Line Mapping

### WSM6 → WDM6 Line Correspondence

| WSM6 Section | Lines | WDM6 Insert After | Content |
|--------------|-------|-------------------|---------|
| Saturation helpers | 72-165 | Line 70 | qvs_water, qvs_ice, fpvs |
| Terminal velocity | 166-250 | After saturation | fall_qi, fall_qs, fall_qg |
| Depositional growth | 1400-1600 | Line 595 | Bergeron, vapor deposition |
| Freezing | 1600-1800 | After deposition | Homogeneous, heterogeneous |
| Melting | 1800-2000 | After freezing | Snow/ice/graupel → rain |
| Ice interactions | 2000-2200 | After melting | Aggregation, riming |
| Ice sedimentation | 2200-2400 | Replace 597-637 | Terminal velocity + fall |

---

## Exact Copy Instructions

### Copy Command Template

```bash
# Extract specific lines from WSM6
sed -n '72,165p' ../WSM6/ERF_AdvanceWSM6.cpp > /tmp/saturation.cpp
sed -n '166,250p' ../WSM6/ERF_AdvanceWSM6.cpp > /tmp/terminal_vel.cpp
sed -n '1400,1600p' ../WSM6/ERF_AdvanceWSM6.cpp > /tmp/deposition.cpp
sed -n '1600,1800p' ../WSM6/ERF_AdvanceWSM6.cpp > /tmp/freezing.cpp
sed -n '1800,2000p' ../WSM6/ERF_AdvanceWSM6.cpp > /tmp/melting.cpp
sed -n '2000,2200p' ../WSM6/ERF_AdvanceWSM6.cpp > /tmp/ice_interactions.cpp
sed -n '2200,2400p' ../WSM6/ERF_AdvanceWSM6.cpp > /tmp/sedimentation.cpp

# Now manually paste these into ERF_AdvanceWDM6.cpp at the locations specified
```

Or use your text editor to copy/paste the ranges.

---

## Before/After File Size

**Current:**
- `ERF_AdvanceWDM6.cpp`: 650 lines

**After ice implementation:**
- `ERF_AdvanceWDM6.cpp`: ~2050 lines (+1400 lines)

**For comparison:**
- `ERF_AdvanceWSM6.cpp`: 2508 lines (similar structure)

---

## Verification Checklist

After implementation, verify:

- [ ] File compiles with no errors
- [ ] File is ~2000-2100 lines (was 650)
- [ ] Search for "Step X" comments - should have ~15-20 steps
- [ ] Search for "qi_arr" - should find 50+ occurrences
- [ ] Search for "qs_arr" - should find 50+ occurrences  
- [ ] Search for "qg_arr" - should find 50+ occurrences
- [ ] Search for "AMREX_GPU_DEVICE" - should find 50+ annotations
- [ ] Search for "ParallelFor" - should find ~25 loops
- [ ] No compilation warnings about unused variables
- [ ] Test run produces qi, qs, qg > 0

---

## Troubleshooting

### Issue: "Cannot find WSM6 ice code at specified lines"

WSM6 file structure may have changed. Search for keywords:
```bash
grep -n "Bergeron" ../WSM6/ERF_AdvanceWSM6.cpp
grep -n "depositional" ../WSM6/ERF_AdvanceWSM6.cpp
grep -n "piacr" ../WSM6/ERF_AdvanceWSM6.cpp  # rain freezing
grep -n "psmlt" ../WSM6/ERF_AdvanceWSM6.cpp  # snow melting
```

### Issue: "Compilation errors after copying"

Check:
1. All helper functions copied (saturation, terminal velocity)
2. No missing braces `}` from copy/paste
3. All array names match (`qi_arr`, `qs_arr`, `qg_arr`)
4. `AMREX_GPU_DEVICE` annotations present

### Issue: "Code compiles but qi/qs/qg stay zero"

Check:
1. Ice nucleation process is included
2. Temperature goes below 0°C in test case
3. Supersaturation w.r.t. ice is positive
4. Depositional growth kernel is being called

---

## Commit After Completion

```bash
git add Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp
git commit -m "[WDM6] Add complete ice physics from WSM6

Ported ice processes from WSM6 C++ GPU implementation:
- Saturation functions (qvs_ice, qvs_water)
- Terminal velocity (fall_qi, fall_qs, fall_qg)
- Depositional growth (Bergeron, vapor deposition)
- Freezing processes (homogeneous, heterogeneous, immersion)
- Melting processes (snow→rain, ice→rain, graupel→rain)
- Ice-ice interactions (aggregation, riming, collection)
- Ice sedimentation (qi, qs, qg terminal velocity + fall)

Implementation: ~1400 lines copied from WSM6
- Ice physics identical (both single-moment)
- Only warm rain differs (WDM6 is double-moment)

File size: 650 → 2050 lines

Status: Complete GPU WDM6 implementation
- Warm rain: double-moment (qc, qr with nc, nr)
- Ice: single-moment (qi, qs, qg)
- All processes on GPU with ParallelFor

Build: cmake -DERF_ENABLE_WDM6_FORT=OFF for full GPU
Test: Run supercell or mountain wave with ice

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

## Quick Reference

**Key file:** `Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp`

**Current line count:** 650

**Target line count:** ~2050

**Main task:** Copy ~1400 lines from WSM6 lines 72-2400

**Est. time:** 2-3 hours

**Difficulty:** Medium (mostly copy/paste)

---

## Questions / Help

If you get stuck, check:
1. This guide - step by step instructions
2. `WDM6_CPP_IMPLEMENTATION_SUMMARY.md` - design decisions
3. `../WSM6/ERF_AdvanceWSM6.cpp` - working reference code
4. Git history - see what was done for warm rain

Good luck! The ice physics implementation is straightforward since it's a direct copy from WSM6.

---

**Last updated:** 2026-07-29 after warm rain GPU implementation
**Status:** Ready to implement ice physics
**Next step:** Follow Step 1 above
