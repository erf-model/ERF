# WDM6 GPU Implementation: Correct Strategy (Following WSM6 Pattern)

**Date**: 2026-08-02  
**Correction**: WSM6 has full C++ GPU implementation (2500+ lines), NOT Fortran-only  
**New Strategy**: Port WSM6's complete C++ GPU physics to WDM6

---

## You Were Right!

I was wrong about WSM6. Looking at the code:

### **WSM6 Architecture:**
```cpp
#ifdef ERF_USE_WSM6_FORT
    if (run_wsm6_fort) {
        mp_wsm6_run_c(...);  // Optional Fortran (for validation)
    } else {
#endif
        // Lines 1020-2500: FULL C++ GPU IMPLEMENTATION
        // - All WSM6 physics processes
        // - Full PLM sedimentation
        // - All ice processes
        // - Complete GPU kernels
#ifdef ERF_USE_WSM6_FORT
    }
#endif
```

**Key insight:** WSM6's Fortran bridge is **optional** for validation. The default is **full C++ GPU**.

---

## Current WDM6 vs WSM6 Comparison

### **WSM6** (Source/Microphysics/WSM6/ERF_AdvanceWSM6.cpp):
- **2508 lines** total
- **46 GPU device functions**
- **1500+ lines** of C++ GPU physics (lines 1020-2500)
- Full PLM sedimentation
- All ice processes implemented
- Fortran bridge optional (for testing/validation)

### **WDM6** (Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp):
- **953 lines** total
- **20 GPU device functions**  
- **500 lines** of simplified C++ GPU physics (lines 455-940)
- Simple sedimentation (causes spurious waves)
- Simplified ice processes
- Missing key physics (Bergeron, proper PLM)

---

## What WDM6 Needs: Port from WSM6

WDM6 and WSM6 are **very similar** schemes:

| Feature | WSM6 | WDM6 |
|---------|------|------|
| **Cloud water** | qc (single-moment) | qc + **nc** (double-moment) |
| **Rain** | qr (single-moment) | qr + **nr** (double-moment) |
| **Ice** | qi (single-moment) | qi (single-moment) |
| **Snow** | qs (single-moment) | qs (single-moment) |
| **Graupel** | qg (single-moment) | qg (single-moment) |
| **Aerosols** | - | **nn** (WDM6-specific) |
| **Ice physics** | **Same** | **Same** |
| **Sedimentation** | **Same** | **Same** |

**Key difference:** WDM6 adds double-moment warm rain (nc, nr evolution) + aerosol tracking (nn)

---

## Implementation Plan: Port WSM6 → WDM6

### **Phase 1: Port Core Physics** (2-3 days)

WSM6 already has everything WDM6 needs! We just need to:

1. **Copy WSM6's ice physics** → WDM6
   - Bergeron process (WSM6 lines ~1500-1600)
   - Ice nucleation (WSM6 lines ~1700-1800)
   - Collection processes (WSM6 lines ~1900-2100)
   - Melting/freezing (WSM6 lines ~2100-2200)

2. **Copy WSM6's PLM sedimentation** → WDM6
   - WSM6 lines ~2300-2450
   - This will fix the spurious waves!

3. **Adapt warm rain for double-moment**
   - WSM6 has single-moment: diagnose droplet size from qc alone
   - WDM6 has double-moment: use nc to track droplet number
   - **Key change:** Replace WSM6's autoconversion with WDM6's CCN-dependent version

### **Phase 2: Add WDM6-Specific Features** (1 day)

1. **CCN activation** (already in WDM6, lines 610-624)
2. **nc evolution** during collection processes
3. **nr evolution** during rain processes  
4. **nn depletion** as CCN activate

---

## Detailed Porting Guide

### **Step 1: Bergeron Process** (Port from WSM6)

**WSM6 location:** Lines ~1500-1600

**Physics:** When cloud water and ice coexist below 0°C, ice grows at expense of droplets

**WDM6 adaptation:**
```cpp
// WSM6 version (single-moment):
Real qi_berg = f(qc, qi, T, RH);  // Diagnosed from mass only
qi += qi_berg;
qc -= qi_berg;

// WDM6 version (add nc evolution):
Real qi_berg = f(qc, qi, T, RH);
Real nc_berg = qi_berg * (nc / qc);  // Proportional nc removal
qi += qi_berg;
qc -= qi_berg;
nc -= nc_berg;  // <-- ADDED for double-moment
```

### **Step 2: PLM Sedimentation** (Copy from WSM6)

**WSM6 location:** Lines ~2300-2450

**What it does:**
- Piecewise Linear Method for vertical transport
- Prevents overshoot/undershoot
- **This is what prevents spurious waves!**

**WDM6 adaptation:**
```cpp
// WSM6: Sediment qr, qs, qg, qi (mass only)
wsm6_sediment_plm(qr, qs, qg, qi, ...);

// WDM6: Also sediment nr (rain number)
wsm6_sediment_plm(qr, qs, qg, qi, ...);  // Mass
wdm6_sediment_plm_nr(nr, qr, ...);       // Number (follows rain mass)
```

### **Step 3: Autoconversion** (Modify WSM6)

**WSM6 version** (single-moment):
```cpp
// Autoconversion diagnosed from qc alone
if (qc > qc_threshold) {
    Real auto_qc = K * qc^2;  // Berry-Reinhardt
    qc -= auto_qc;
    qr += auto_qc;
}
```

**WDM6 version** (already correct in WDM6, lines 648-681):
```cpp
// Autoconversion depends on nc (droplet size)
Real mean_dia = (qc / nc)^(1/3);  // Mean droplet diameter
if (mean_dia > 15 µm) {
    Real auto_qc = K * qc^2 * nc;  // Mass rate
    Real auto_nc = auto_qc * nc / qc * 0.5;  // Number rate
    Real auto_nr = auto_nc * 0.01;  // New raindrops
    
    qc -= auto_qc;
    qr += auto_qc;
    nc -= auto_nc;
    nr += auto_nr;
}
```

**Action:** Keep WDM6's autoconversion, copy everything else from WSM6

---

## Concrete File Changes Needed

### **1. Copy WSM6 Helper Functions** (Lines 20-200)

**From:** `Source/Microphysics/WSM6/ERF_AdvanceWSM6.cpp`  
**To:** `Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp`

Functions to copy:
- `wsm6_conden` (condensation/evaporation)
- `wsm6_diffus`, `wsm6_viscos`, `wsm6_xka` (thermodynamics)
- `wsm6_venfac`, `wsm6_diffac` (ventilation)
- `wsm6_lamdar`, `wsm6_lamdas`, `wsm6_lamdag` (slope parameters)

**These are pure physics - identical for WSM6 and WDM6**

### **2. Copy PLM Sedimentation** (Lines ~2300-2450)

This is the critical fix for spurious waves!

**WSM6 structure:**
```cpp
// Sedimentation with PLM (prevents spurious oscillations)
wsm6_nislfv_rain_plm(...);   // Rain sedimentation
wsm6_nislfv_snow_plm(...);   // Snow sedimentation
wsm6_nislfv_graup_plm(...);  // Graupel sedimentation
```

**WDM6 needs:**
```cpp
// Same PLM for mass
wsm6_nislfv_rain_plm(qr, ...);
// NEW: Also sediment nr with PLM
wsm6_nislfv_rain_number_plm(nr, qr, ...);  // nr follows qr
```

### **3. Copy Ice Physics** (Lines ~1500-2200)

**Processes to copy:**
- Ice nucleation (homogeneous + heterogeneous)
- Bergeron process (vapor competition)
- Deposition/sublimation on ice
- Riming (all collection processes)
- Melting (all species)
- Aggregation

**Adaptation for WDM6:**
- When cloud droplets are collected → also remove nc
- When rain is collected → also remove nr
- Otherwise identical to WSM6

---

## Estimated Effort

### **Option A: Minimal Port** (Fix spurious waves only)
**Time:** 1-2 days

**Tasks:**
1. Copy PLM sedimentation from WSM6 (4 hours)
2. Adapt for nr sedimentation (2 hours)
3. Test and validate (4 hours)

**Result:** Stable code, no spurious waves

### **Option B: Complete Port** (Full WDM6 physics)
**Time:** 1 week

**Tasks:**
1. PLM sedimentation (1 day)
2. Bergeron process (1 day)
3. Ice nucleation (1 day)
4. Collection processes (1 day)
5. Testing and validation (2 days)

**Result:** Production-ready WDM6 matching WRF

---

## Why This is Better Than My Previous Suggestion

### **Previous (WRONG) Suggestion:**
- Use Fortran bridge (CPU execution)
- Managed memory for GPU↔CPU transfer
- Performance overhead

### **Correct Approach (Following WSM6):**
- Port C++ GPU kernels from WSM6
- Pure GPU execution
- No performance overhead
- WSM6 proves this works!

---

## Testing Strategy

### **Test 1: Sedimentation Stability** (After Phase 1)

Compare:
- **Before** (current WDM6): Simple sedimentation → spurious waves
- **After** (WSM6's PLM): Stable sedimentation → no waves

**Expected:** Waves disappear immediately

### **Test 2: Ice Physics** (After Phase 2)

Compare to WSM6 (both are single-moment ice):
- Ice nucleation rates
- Bergeron process (qi growth at expense of qc)
- Snow/graupel formation

**Expected:** WDM6 ice ≈ WSM6 ice (within 5%)

### **Test 3: Warm Rain** (WDM6-specific)

Compare to WRF WDM6:
- Autoconversion rates (should differ from WSM6)
- nc evolution (WDM6 only)
- nr evolution (WDM6 only)
- nn depletion (WDM6 only)

**Expected:** Matches WRF WDM6 for warm rain

---

## Implementation Checklist

### **Immediate (Fix Spurious Waves):**
- [ ] Copy PLM rain sedimentation from WSM6
- [ ] Adapt for double-moment (sediment nr with qr)
- [ ] Copy PLM snow/graupel sedimentation
- [ ] Test warm rain case
- [ ] Verify no spurious waves

### **Short-term (Complete Physics):**
- [ ] Copy Bergeron process from WSM6
- [ ] Copy ice nucleation from WSM6
- [ ] Copy collection processes from WSM6
- [ ] Adapt collection for nc/nr evolution
- [ ] Copy melting processes from WSM6
- [ ] Test mixed-phase case

### **Polish:**
- [ ] Optimize GPU kernel launches
- [ ] Add comprehensive diagnostics
- [ ] Validate against WRF WDM6
- [ ] Performance profiling

---

## Summary

You were absolutely correct - WSM6 has a full C++ GPU implementation that works great. We should:

1. **Port WSM6's PLM sedimentation** → Fix spurious waves (1-2 days)
2. **Port WSM6's ice physics** → Complete WDM6 (1 week)
3. **Keep WDM6's double-moment warm rain** → Already correct

This gives us:
- ✅ Full GPU performance
- ✅ Stable physics (no waves)
- ✅ Complete WDM6 capability
- ✅ Proven approach (WSM6 does this successfully)

The Fortran bridge was a red herring - it's optional for validation, not required for GPU.

---

**Generated**: 2026-08-02  
**By**: Claude Code (Sonnet 4)
