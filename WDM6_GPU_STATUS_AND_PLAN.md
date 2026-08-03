# WDM6 GPU Implementation: Status and Action Plan

**Date**: 2026-08-02  
**Status**: Warm rain + simplified ice implemented, spurious waves reported  
**Goal**: Complete full WDM6 GPU implementation for production use

---

## Current Implementation Status

### ✅ **IMPLEMENTED (C++ GPU Kernels)**

#### **Warm Rain Physics (Double-Moment):**
- ✅ CCN activation (WDM6-specific, lines 610-624)
- ✅ Condensation/evaporation (lines 629-643)
- ✅ Autoconversion qc→qr with nc evolution (lines 646-681)
- ✅ Accretion (rain collecting cloud water, lines 683-707)
- ✅ Rain evaporation with nr evolution (lines 709-728)
- ✅ Rain sedimentation (double-moment with nr, lines 853-872)

#### **Ice Physics (Simplified Single-Moment):**
- ✅ Homogeneous freezing (<-40°C, line 740-748)
- ✅ Heterogeneous freezing (0 to -40°C, lines 750-758)
- ✅ Rain freezing → graupel (lines 760-770)
- ✅ Ice → snow aggregation (lines 772-778)
- ✅ Snow melting → rain (lines 780-788)
- ✅ Graupel melting → rain (lines 790-797)
- ✅ Riming (snow + cloud → graupel, lines 799-809)
- ✅ Snow sedimentation (single-moment, lines 874-894)
- ✅ Graupel sedimentation (single-moment, lines 896-910)
- ✅ Ice sedimentation (simplified, lines 912-924)

#### **Infrastructure:**
- ✅ Saturation calculations (water and ice, lines 559-591)
- ✅ Terminal velocity helpers (lines 24-196)
- ✅ Density factor corrections (lines 554-556)
- ✅ Theta↔temperature conversion (CRITICAL BUG FIX, lines 439-453)
- ✅ Minor timestep sub-cycling (lines 550-826)
- ✅ Global diagnostics (precipitation, mixing ratios)

---

## ❌ **NOT YET IMPLEMENTED (Missing Physics)**

### **Major WRF WDM6 Processes Not in C++ GPU Code:**

#### **1. Bergeron Process (Vapor Deposition on Ice)**
**What it does:** When both cloud water and ice coexist below 0°C, ice grows at expense of droplets via vapor pathway (ice has lower saturation vapor pressure)

**WRF Location:** `module_mp_wdm6.F` lines ~2500-2600  
**Impact:** Critical for mixed-phase clouds. Without it:
- Ice crystals grow too slowly
- Cloud water persists unrealistically at cold temperatures
- Snow/precipitation underestimated

**Priority:** 🔴 **HIGH**

---

#### **2. Ice Number Concentration Tracking (ni)**
**What it does:** WRF WDM6 *optionally* tracks ice crystal number concentration for better size distribution

**Current C++ approach:** Single-moment ice (fixed N₀)  
**WRF approach:** Can use double-moment ice when `QNIN > 0` in wrfinput

**Impact:** Lower fidelity ice microphysics  
**Priority:** 🟡 **MEDIUM** (single-moment is acceptable for many applications)

---

#### **3. Detailed Terminal Velocity Calculations**
**What's simplified:**
- Rain: Using simplified lambda calculation (line 855-859)
- Snow: Temperature-dependent N₀ but simplified (lines 876-884)
- Graupel: Simplified fixed N₀ (lines 898-900)
- Ice: Constant fall speed (line 915)

**WRF approach:** Full slope parameter calculations with diameter limits, temperature dependencies, density corrections

**Impact:** 
- Sedimentation rates may differ by 10-30%
- Affects precipitation timing and amount
- Particle size distributions less accurate

**Priority:** 🟡 **MEDIUM**

---

#### **4. Piecewise Linear Method (PLM) Sedimentation**
**Current C++ approach:** Top-down simple flux (lines 843-931)

**WRF approach:** Full PLM with limiter to prevent overshoots

**Impact:** 
- Can cause spurious oscillations (THIS MAY BE YOUR SPURIOUS WAVES!)
- Numerical diffusion in precipitation
- Less accurate vertical transport

**Priority:** 🔴 **HIGH** (likely cause of spurious waves)

---

#### **5. Graupel/Snow Collection Processes**
**Missing:**
- Snow collecting rain → graupel
- Graupel collecting cloud water (riming)
- Snow self-collection (aggregation rate)
- Graupel collecting snow

**Current:** Only have snow+cloud→graupel riming (line 800-809)

**Impact:** Ice growth rates underestimated  
**Priority:** 🟡 **MEDIUM**

---

#### **6. Ice Multiplication (Hallett-Mossop Process)**
**What it does:** Secondary ice production via rime splintering at -3 to -8°C

**WRF Location:** `module_mp_wdm6.F` lines ~2700-2800  
**Impact:** Ice number concentration enhancement in specific T range  
**Priority:** 🟢 **LOW** (important for specific cloud regimes, but not fundamental)

---

#### **7. Saturation Adjustment Over Ice**
**Current:** Only doing condensation/evaporation for T > 0°C (line 630)

**Missing:** Vapor deposition/sublimation for ice saturation at T < 0°C

**Impact:** 
- Ice crystals don't grow from supersaturated vapor
- Incorrectly handles mixed-phase conditions
- Related to missing Bergeron process

**Priority:** 🔴 **HIGH**

---

## 🐛 **Reported Issues**

### **1. Spurious Waves with Warm Rain Only**
**User Report:** "lot of unphysical spurious waves"

**Likely Causes:**
1. **Simple sedimentation scheme** (no PLM) causing numerical oscillations
2. **Missing limiter** on fall speeds / fluxes
3. **Latent heating feedback** to dynamics (theta conversion is correct now, but rate might be too large per timestep)
4. **Minor timestep too large** (currently hardcoded to 2, line 548)

**Debugging Steps:**
- Reduce minor timesteps (try 4-6 instead of 2)
- Check CFL condition on sedimentation: `v_term * dt / dz < 1`
- Add limiters on theta updates per timestep
- Implement PLM sedimentation

---

## 📋 **Implementation Priority Roadmap**

### **Phase 1: Fix Spurious Waves (IMMEDIATE)**

**Goal:** Make current warm-rain code stable and production-ready

**Tasks:**
1. ✅ Verify theta conversion is working (ALREADY DONE)
2. ⬜ Implement PLM sedimentation for rain
   - Location: Replace lines 843-931
   - Reference: ERF Morrison has PLM implementation
   - Prevents overshoot/undershoot oscillations
3. ⬜ Add CFL limiter on sedimentation velocities
4. ⬜ Increase minor timestep count (line 548: `wdm6_loops = max(2, ceil(dt/30))`)
5. ⬜ Add latent heat rate limiter (cap theta change per minor timestep)

**Expected Result:** Stable warm-rain simulations without spurious waves

---

### **Phase 2: Complete Ice Physics (NEXT)**

**Goal:** Full WDM6 capability matching WRF

**Tasks:**
1. ⬜ **Bergeron process** (vapor competition between water and ice)
   - Critical for mixed-phase clouds
   - Lines ~2500-2600 in WRF `module_mp_wdm6.F`
   - Implement as Step 8.5 between current Steps 8 and 9
   
2. ⬜ **Saturation adjustment for ice** (T < 0°C)
   - Modify Step 5 (lines 628-643) to handle ice deposition/sublimation
   - Use `qsati_arr` (already computed) for ice saturation
   
3. ⬜ **Ice-snow-graupel collection processes**
   - Snow collecting rain
   - Snow self-collection
   - Graupel collecting snow
   - Reference: WRF lines ~2650-2750

**Expected Result:** Realistic mixed-phase cloud behavior, improved precipitation

---

### **Phase 3: Terminal Velocity Refinement (POLISH)**

**Goal:** Match WRF sedimentation rates

**Tasks:**
1. ⬜ Implement full slope parameter calculations (lambdar, lamdas, lamdag)
2. ⬜ Add diameter limiters (dicon, dimin, dimax)
3. ⬜ Temperature-dependent N₀ for snow (already partially done, refine)
4. ⬜ Density-weighted fall speeds for all species

**Expected Result:** Precipitation timing and amounts match WRF within 5%

---

### **Phase 4: Optional Enhancements (FUTURE)**

1. ⬜ Double-moment ice (track ni)
2. ⬜ Ice multiplication (Hallett-Mossop)
3. ⬜ Advanced CCN activation with vertical velocity coupling

---

## 🔍 **Immediate Action: Diagnose Spurious Waves**

Before implementing new physics, let's understand the current issue:

### **Diagnostic Steps:**

1. **Check sedimentation CFL:**
   ```cpp
   // In sedimentation loop, add diagnostic:
   Real cfl = vt * dt_advance / dz_sedi;
   if (cfl > 1.0) {
       amrex::Print() << "CFL violation: cfl=" << cfl 
                     << " at i,j,k=" << i << "," << j << "," << kk << "\n";
   }
   ```

2. **Check theta changes:**
   ```cpp
   // Before/after microphysics:
   Real theta_change = theta_arr(i,j,k) - theta_before;
   if (abs(theta_change) > 1.0) {  // More than 1K change per step
       amrex::Print() << "Large theta change: dtheta=" << theta_change << " K\n";
   }
   ```

3. **Visualize vertical velocity after microphysics:**
   - Spurious waves often show up as alternating updraft/downdraft patterns
   - Check w-velocity field correlation with precipitation regions

4. **Test with smaller timestep:**
   - If waves disappear with dt/2, it's a CFL/stability issue
   - If waves persist, it's likely a physics formulation issue

---

## 🛠 **How to Build GPU Code**

### **Current Build:**
```bash
cd /g/g10/wise14/compiling/clean/ERF
cmake -B build -DERF_ENABLE_WDM6_CPP=ON -DAMREX_GPU_BACKEND=CUDA
cmake --build build -j8
```

### **Test Fortran Comparison:**
```bash
cmake -B build -DERF_ENABLE_WDM6_FORT=ON  # Uses Fortran version
```

The code automatically selects:
- `ERF_USE_WDM6_FORT` defined → Fortran (CPU, full physics)
- Otherwise → C++ GPU kernels (current implementation)

---

## 📊 **Testing Strategy**

### **1. Warm Rain Test Case (Verify Fix):**
- Temperature > 0°C everywhere
- Only qc, qr, nc, nr, nn active
- Should be stable with no spurious waves after Phase 1

### **2. Mixed-Phase Test Case (After Phase 2):**
- Temperature -20 to +10°C
- All species active (qi, qs, qg)
- Compare to WRF: cloud structure, precipitation partition

### **3. Deep Convection Test Case (Full Test):**
- After Phase 3 complete
- Strong updrafts, full ice processes
- Compare surface rain/snow accumulation

---

## 📁 **Key Files**

| File | Purpose | Lines to Focus |
|------|---------|----------------|
| `ERF_AdvanceWDM6.cpp` | Main GPU kernels | 550-940 (physics loop) |
| `ERF_InitWDM6.cpp` | Initialization | 39-189 (nn init, Copy_State) |
| `ERF_UpdateWDM6.cpp` | State synchronization | Full file |
| `ERF_WDM6.H` | Class definition | 19-38 (MicVar enum) |
| WRF `module_mp_wdm6.F` | Reference implementation | 2300-3500 (wdm62D) |

---

## 🎯 **Immediate Next Steps**

1. **Diagnose spurious waves** (add CFL/theta diagnostics)
2. **Implement PLM sedimentation** for rain (lines 843-931)
3. **Test warm-rain stability** (should eliminate waves)
4. **Implement Bergeron process** (critical for ice growth)
5. **Add ice saturation adjustment**
6. **Full mixed-phase testing**

**Estimated Time:**
- Phase 1 (fix waves): 1-2 days
- Phase 2 (complete ice): 3-4 days
- Phase 3 (refinement): 2-3 days
- **Total: ~1-2 weeks for production-ready GPU WDM6**

---

**Generated**: 2026-08-02  
**By**: Claude Code (Sonnet 4)
