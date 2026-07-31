# WDM6 Default Parameter Issues - Why ERF ≠ WRF

**Date**: 2026-07-31  
**Issue**: ERF WDM6 produces different results than WRF WDM6 with "default" settings  
**Root Cause**: ERF has several non-WRF-default behaviors

---

## **Critical Issues Found**

### **Issue #1: Land Mask Hardcoded to "Land"** ⚠️ HIGH PRIORITY

**ERF** (`ERF_AdvanceWDM6.cpp:380`):
```cpp
xland_arr(i,j,k) = Real(1.0);  // ALWAYS LAND!
```

**WRF**:
```fortran
xland(i,j) from WPS land/water mask
  = 1.0 over land
  = 2.0 over water
```

**Impact**:
- Maritime clouds treated as continental
- Wrong CCN activation parameters
- Different autoconversion rates

**Fix Required**: Connect to ERF's `lmask_lev` or initialize from input

---

### **Issue #2: NC Initialization Forced** ⚠️ MEDIUM PRIORITY

**ERF** (`ERF_InitWDM6.cpp:144-153`):
```cpp
if (nc(i,j,k) < 1e1 || nc_physical > 5e8) {
    if (qc(i,j,k) > 1e-9) {
        nc(i,j,k) = 5.0e7 / rho(i,j,k);  // Forces 50 cm⁻³
    } else {
        nc(i,j,k) = 1.e1;  // ncmin
    }
}
```

**WRF**:
```fortran
! NO forced initialization of nc
! nc starts at zero or from restart
! Built up naturally by CCN activation
```

**Impact**:
- First timestep has artificial nc=50 cm⁻³ where qc > 0
- Skips natural CCN activation spin-up
- May cause transient artifacts

**Fix Required**: Remove forced initialization, let physics build nc naturally

---

### **Issue #3: NN Re-initialization Every Timestep** ⚠️ MEDIUM PRIORITY

**ERF** (`ERF_InitWDM6.cpp:130-132`):
```cpp
// Called EVERY Copy_State_to_Micro (every timestep)
if (nn(i,j,k) < Real(1.0)) {
    nn(i,j,k) = ccn0_local / rho(i,j,k);
}
```

**WRF** (`module_mp_wdm6.F:219-227`):
```fortran
! Called ONLY on itimestep == 1
if (itimestep .eq. 1) then 
  do j = jms,jme
    do k = kms,kme    
      do i = ims,ime
        nn(i,k,j) = ccn0   
      enddo
    enddo
  enddo
endif
```

**Impact**:
- ERF replenishes nn every timestep if it drops below 1.0 #/kg
- WRF initializes once, then nn evolves freely
- Different aerosol depletion behavior

**Fix Required**: Initialize nn only on first call, not every timestep

---

### **Issue #4: CCN0 Units** ⚠️ NEEDS VERIFICATION

**Question**: Does WRF pass ccn0 in #/m³ or #/kg?

**ERF Assumption**:
```cpp
Real m_ccn0{100.0e6};  // #/m³
nn = m_ccn0 / rho;     // Convert to #/kg
```

**WRF Code**:
```fortran
nn(i,k,j) = ccn0  ! Direct assignment
```

**Possible scenarios**:
1. WRF passes ccn0 in #/m³, ERF should divide by ρ ✅ (ERF correct)
2. WRF passes ccn0 in #/kg, ERF should NOT divide by ρ ❌ (ERF wrong)

**Status**: Need to check WRF namelist and driver to confirm units

---

### **Issue #5: NCmin Threshold for Re-init** ⚠️ LOW PRIORITY

**ERF**: Checks `nc < 1e1` OR `nc_physical > 5e8` for re-initialization

**WRF**: No such check, nc can be any value

**Impact**: Minor, but causes ERF to "fix" low nc values WRF wouldn't touch

---

## **Recommended Fixes**

### **Fix #1: Remove Land Mask Hardcoding** (Immediate)

```cpp
// In ERF_AdvanceWDM6.cpp line 374-381
// OLD:
FArrayBox xland_fab(Box(IntVect(imlo,jmlo,0), IntVect(imhi,jmhi,0)), 1, The_Pinned_Arena());
auto const& xland_arr = xland_fab.array();
ParallelFor(Box(IntVect(imlo,jmlo,0), IntVect(imhi,jmhi,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    xland_arr(i,j,k) = Real(1.0);  // Default to land
});

// NEW:
// TODO: Get actual land mask from ERF's m_lmask
// For now, read from input file or use ERF::lmask_lev
// If not available, default to water (2.0) for maritime conditions
FArrayBox xland_fab(Box(IntVect(imlo,jmlo,0), IntVect(imhi,jmhi,0)), 1, The_Pinned_Arena());
auto const& xland_arr = xland_fab.array();
ParallelFor(Box(IntVect(imlo,jmlo,0), IntVect(imhi,jmhi,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    xland_arr(i,j,k) = Real(2.0);  // Default to WATER (maritime)
});
```

**Rationale**: Maritime (water) is more common for idealized cases. If user wants continental, they should set it explicitly.

### **Fix #2: Remove Forced NC Initialization** (Important)

```cpp
// In ERF_InitWDM6.cpp line 134-153
// DELETE THIS ENTIRE BLOCK:
// Real nc_physical = nc(i,j,k) * rho(i,j,k);
// Real nc_max_physical = Real(5.0e8);
// if (nc(i,j,k) < Real(1.e1) || nc_physical > nc_max_physical) {
//     if (qc(i,j,k) > Real(1.e-9)) {
//         nc(i,j,k) = Real(5.0e7) / rho(i,j,k);
//     } else {
//         nc(i,j,k) = Real(1.e1);
//     }
// }

// REPLACE WITH:
// Only enforce absolute minimum (like WRF)
nc(i,j,k) = amrex::max(nc(i,j,k), Real(1.e1));  // WRF's ncmin
```

**Rationale**: Let CCN activation build nc naturally from zero, like WRF does.

### **Fix #3: Initialize NN Only Once** (Important)

```cpp
// In ERF_InitWDM6.cpp, add static flag:
static bool nn_initialized = false;

// In Copy_State_to_Micro, around line 130:
if (!nn_initialized) {
    ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (nn(i,j,k) < Real(1.0)) {
            nn(i,j,k) = ccn0_local / rho(i,j,k);
        }
    });
    nn_initialized = true;
}
```

**Rationale**: Match WRF's `itimestep == 1` behavior - initialize once, then let nn evolve.

### **Fix #4: Verify CCN0 Units** (Investigation)

Check WRF driver code to confirm:
1. What units does WRF pass for ccn0?
2. Does WRF's initialization do `nn = ccn0` or `nn = ccn0/rho`?

If WRF passes #/kg directly, then ERF should do:
```cpp
nn(i,j,k) = ccn0_local;  // No division by rho
```

---

## **Testing After Fixes**

### **Test Case**: WRF Validation

Run identical case in WRF and ERF:
1. Same initial conditions (T, P, qv, qc=0, qi=0, etc.)
2. Same domain size and resolution
3. Same timestep
4. Compare after 100 timesteps:
   - nc, nn, nr fields
   - Surface precipitation
   - Cloud/rain distribution

### **Expected After Fixes**:

✅ nc starts at zero, builds up via CCN activation  
✅ nn initializes to ccn0/ρ once, then depletes naturally  
✅ xland uses actual land/water mask (or maritime default)  
✅ Results match WRF to within numerical roundoff

---

## **Summary**

ERF WDM6 has **4 significant deviations** from WRF defaults:

| Issue | Priority | Impact | Status |
|-------|----------|--------|--------|
| Land mask = 1.0 (land) | HIGH | Wrong CCN activation | Easy fix |
| Forced nc = 50 cm⁻³ | MEDIUM | Skips natural spin-up | Easy fix |
| NN re-init every step | MEDIUM | Different aerosol evolution | Medium fix |
| CCN0 units unclear | TBD | May be critical | Needs investigation |

**Recommendation**: Fix #1, #2, #3 immediately. Investigate #4 before finalizing.

---

**Generated**: 2026-07-31  
**By**: Claude Code (Sonnet 4)
