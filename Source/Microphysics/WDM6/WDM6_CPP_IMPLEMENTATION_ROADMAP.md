# WDM6 C++ GPU Implementation Roadmap

**Strategy:** Adapt existing WSM6 C++ GPU code instead of porting WRF Fortran from scratch

**Cost:** ~550 new lines vs 3237 lines (5x cheaper)

**Timeline:** 4-6 hours with testing, 2-3 hours for core implementation

**Reference Files:**
- Source: `/g/g10/wise14/compiling/clean/ERF/Source/Microphysics/WSM6/ERF_AdvanceWSM6.cpp` (2508 lines)
- Target: `/g/g10/wise14/compiling/clean/ERF/Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp` (currently 257 lines)
- WRF Reference: `/p/lustre1/wise14/wrf/compiling/WRF/WRF_v4.7.1_clean/phys/module_mp_wdm6.F` (3237 lines)

---

## Phase Overview

### Phase 2A: Infrastructure (30-45 min, ~20k tokens)
Add number concentration arrays and basic structure

### Phase 2B: CCN Activation (45-60 min, ~25k tokens)
Port WRF activ_conc to C++ GPU device function

### Phase 2C: Autoconversion (30-45 min, ~15k tokens)
Modify cloud→rain conversion for double-moment

### Phase 2D: Sedimentation (30-45 min, ~15k tokens)
Couple mass and number concentration tracking

### Phase 2E: Testing & Debug (1-2 hours, ~20k tokens)
Compile, fix errors, validate against WSM6

**Total Estimated:** 4-6 hours, ~95k tokens

---

## What WSM6 Code Can Be Reused (80%)

### ✅ No Changes Needed - Copy Directly

1. **Ice processes** (lines 1146-2400 in WSM6)
   - Bergeron process
   - Ice crystal growth
   - Freezing/melting
   - Graupel formation
   - **Action:** Copy as-is

2. **Saturation calculations** (lines 1146-1200)
   - qvs, qvsi functions
   - Relative humidity
   - **Action:** Copy as-is

3. **Terminal velocity functions** (lines 1201-1300)
   - Fall speeds for rain, snow, graupel
   - **Action:** Copy as-is (WDM6 uses same formulas)

4. **GPU infrastructure**
   - ParallelFor loops
   - Device function patterns
   - FArrayBox allocations
   - **Action:** Copy structure

5. **Helper functions** (already in ERF_AdvanceWDM6.cpp lines 20-67)
   - wdm6_cpmcal, wdm6_xlcal, etc.
   - **Action:** Already done ✅

---

## What Needs to Change (20%)

### ⚠️ Add Number Concentration Arrays

**Where:** Throughout Advance() function

**What to add:**

```cpp
// After line 70 in ERF_AdvanceWDM6.cpp, in Advance() function:
auto const& nn_arr = mic_fab_vars[MicVar_WDM6::nn]->array(mfi);  // Aerosol number
auto const& nc_arr = mic_fab_vars[MicVar_WDM6::nc]->array(mfi);  // Cloud droplet number
auto const& nr_arr = mic_fab_vars[MicVar_WDM6::nr]->array(mfi);  // Rain drop number
```

**Impact:** ~50 locations need nn, nc, nr added to function signatures and loops

---

## Detailed Implementation Steps

---

## Step 1: Copy WSM6 Structure (15 min)

### 1.1: Copy Main Loop Structure

From `ERF_AdvanceWSM6.cpp` lines 890-960, copy:
- Physical constants setup
- MFIter loop initialization
- Array pointer extraction
- Box definitions

**Target location:** Replace current stub in `ERF_AdvanceWDM6.cpp` Advance() after line 100

**Changes needed:**
- Change `MicVar_WSM6::` → `MicVar_WDM6::`
- Add nn_arr, nc_arr, nr_arr pointers
- Keep ccn0 parameter (already have)

### 1.2: Copy Working FArrayBox Allocations

From `ERF_AdvanceWSM6.cpp` lines 1025-1095, copy all FArrayBox declarations:

```cpp
FArrayBox denfac_fab(fab_box,1, Arena_Used);
FArrayBox xni_fab(fab_box,1, Arena_Used);
// ... all other working arrays
```

**Target:** After array pointer setup in WDM6 Advance()

**Changes needed:**
- Add 3 new FArrayBaxes for nn, nc, nr if not already in mic_fab_vars
- Otherwise, none

### 1.3: Copy Initialization Loops

From `ERF_AdvanceWSM6.cpp` lines 1146-1200, copy:
- Saturation calculation loop
- Density factor calculation
- Initial state setup

**Target:** After FArrayBox allocations

**Changes needed:**
- Add nc, nr initialization to physical minimums (ncmin, nrmin)
- Initialize nn to ccn0 if not already set

---

## Step 2: Add Number Concentration Tracking (30 min)

### 2.1: Slope Parameter Functions

**WSM6 rain slope (single-moment):**
```cpp
// WSM6 line ~1400
lamdar = cbrt(pidn0r / (den * qr))  // Fixed N0
```

**WDM6 rain slope (double-moment):**
```cpp
// Use actual nr instead of fixed N0
lamdar = cbrt((pidnr * nr) / (den * qr))
```

**Action:** Update slope calculation device functions

**Reference:** WRF module_mp_wdm6.F line 2240:
```fortran
lamdar(x,y,z)= exp(log(((pidnr*z)/(x*y)))*((.33333333)))
```

**Where:** Create new device function around line 70:

```cpp
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_lamdar(Real qr, Real den, Real nr) {
    // qr = rain mixing ratio (kg/kg)
    // den = air density (kg/m^3)
    // nr = rain number concentration (#/m^3)
    return std::cbrt((pidnr * nr) / (den * qr));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_lamdac(Real qc, Real den, Real nc) {
    // qc = cloud mixing ratio (kg/kg)
    // den = air density (kg/m^3)
    // nc = cloud droplet number concentration (#/m^3)
    return std::cbrt((pidnc * nc) / (den * qc));
}
```

### 2.2: Minimum Enforcement

**Add after saturation loop:**

```cpp
ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    // Enforce minimum number concentrations
    nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k), ncmin);
    nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k), nrmin);
    nn_arr(i,j,k) = amrex::max(nn_arr(i,j,k), Real(0.0));
});
```

**Reference:** WRF module_mp_wdm6.F lines 86-87:
```fortran
real, parameter, private :: ncmin = 1.e1
real, parameter, private :: nrmin = 1.e-2
```

---

## Step 3: Add CCN Activation (45 min)

### 3.1: Create CCN Activation Device Function

**Reference:** WRF module_mp_wdm6.F lines 2040-2200 (activ_conc subroutine)

**Location:** Add before Advance() function, around line 68

**Signature:**
```cpp
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wdm6_ccn_activation(
    Real& nc,           // cloud droplet number (in/out)
    Real& nn,           // aerosol number (in/out)
    const Real qv,      // water vapor mixing ratio
    const Real qc,      // cloud water mixing ratio
    const Real temp,    // temperature (K)
    const Real pres,    // pressure (Pa)
    const Real den,     // air density (kg/m^3)
    const Real w,       // vertical velocity (m/s)
    const Real dt,      // timestep (s)
    const Real ccn0,    // background CCN (#/m^3)
    const bool is_water // true=maritime, false=continental
) {
    // Abdul-Razzak & Ghan (2000) CCN activation
    
    constexpr Real satmax = 1.0048;  // Max supersaturation for activation
    constexpr Real actk = 0.6;       // Activation parameter
    constexpr Real actr = 1.5e-6;    // Activated droplet radius (m)
    
    // Calculate supersaturation
    Real qvs = wdm6_qvs(temp, pres, den);  // Saturation mixing ratio
    Real supersaturation = qv / qvs - 1.0;
    
    if (supersaturation > 0.0 && qc > 1.e-6) {
        // Limit supersaturation for numerical stability
        supersaturation = amrex::min(supersaturation, satmax - 1.0);
        
        // Activation: convert aerosols (nn) to cloud droplets (nc)
        // Based on supersaturation and vertical velocity
        Real s_crit = actk * std::pow(w, 0.5);  // Critical supersaturation
        
        if (supersaturation > s_crit) {
            Real nc_activated = nn * (supersaturation - s_crit) / satmax;
            nc_activated = amrex::min(nc_activated, nn);  // Can't activate more than available
            
            // Update concentrations
            nc += nc_activated;
            nn -= nc_activated;
            
            // Enforce bounds
            nc = amrex::max(nc, ncmin);
            nn = amrex::max(nn, Real(0.0));
        }
    }
}
```

**Simplification Note:** This is a **simplified** activation. Full WRF version includes:
- Size distribution integration
- Hygroscopicity parameters
- Multiple aerosol modes
- See WRF lines 2040-2200 for complete implementation

**For initial implementation:** Use simplified version above, enhance later if needed

### 3.2: Call CCN Activation

**Add in main physics loop, after condensation/evaporation:**

```cpp
ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    // ... condensation/evaporation code from WSM6 ...
    
    // WDM6: CCN activation
    Real w_velocity = /* extract from velocity field */;
    bool is_maritime = (xland_arr(i,j,0) == 0);  // 0=water, 1=land
    
    wdm6_ccn_activation(
        nc_arr(i,j,k), nn_arr(i,j,k),
        qv_arr(i,j,k), qc_arr(i,j,k),
        t_arr(i,j,k), p_arr(i,j,k), den_arr(i,j,k),
        w_velocity, dt, ccn0, is_maritime
    );
});
```

**TODO:** Determine how to get vertical velocity `w` from ERF state

---

## Step 4: Modify Autoconversion (30 min)

### 4.1: WSM6 Autoconversion (Single-Moment)

**From ERF_AdvanceWSM6.cpp** (approximate line 1600):
```cpp
// WSM6: qc → qr conversion based only on qc
if (qc > qc_threshold) {
    Real auto_rate = k_auto * qc * qc;  // Kessler-type
    qc -= auto_rate * dt;
    qr += auto_rate * dt;
}
```

### 4.2: WDM6 Autoconversion (Double-Moment)

**Modify to also convert nc → nr:**

```cpp
// WDM6: Convert both mass (qc→qr) and number (nc→nr)
if (qc > qc_threshold && nc > ncmin) {
    // Mass autoconversion (similar to WSM6)
    Real mean_diameter = wdm6_mean_droplet_diameter(qc, nc, den);
    
    if (mean_diameter > di15) {  // 15 micron threshold
        // Autoconversion rate depends on nc and qc
        Real auto_qc = qck1 * qc * qc * nc;  // kg/kg/s
        auto_qc = amrex::min(auto_qc * dt, qc);
        
        // Number conversion: proportional to mass, but with different efficiency
        Real auto_nc = auto_qc * nc / qc * 0.5;  // 50% efficiency
        auto_nc = amrex::min(auto_nc, nc);
        
        // Update mixing ratios and numbers
        qc -= auto_qc;
        qr += auto_qc;
        nc -= auto_nc;
        nr += auto_nc;
        
        // Enforce minimums
        nc = amrex::max(nc, ncmin);
        nr = amrex::max(nr, nrmin);
    }
}
```

**Reference:** WRF module_mp_wdm6.F lines 900-1000 (autoconversion section in wdm62d)

**Key WDM6 Parameters:**
```cpp
constexpr Real qck1 = 4706.08203;  // From WRF line 2115
constexpr Real di15 = 15.e-6;      // 15 micron threshold
```

**Add helper function:**
```cpp
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_mean_droplet_diameter(Real qc, Real nc, Real den) {
    // qc in kg/kg, nc in #/m^3, den in kg/m^3
    if (nc < ncmin || qc < 1.e-9) return 0.0;
    
    constexpr Real rho_water = 1000.0;  // kg/m^3
    constexpr Real pi = 3.14159265359;
    
    // Volume-weighted mean diameter
    Real volume = qc * den / (rho_water * nc);
    Real diameter = std::pow(6.0 * volume / pi, 1.0/3.0);
    
    return diameter;
}
```

---

## Step 5: Update Sedimentation (45 min)

### 5.1: WSM6 Sedimentation (Single-Moment)

**From ERF_AdvanceWSM6.cpp** (lines 2000-2400):
- Rain falls at terminal velocity based on qr only
- Updates qr at each level
- Accumulates surface precipitation

### 5.2: WDM6 Sedimentation (Double-Moment)

**Key difference:** Must sediment nr along with qr, maintaining their ratio

**Strategy:**
1. Calculate terminal velocity from qr AND nr (not just qr)
2. Sediment qr and nr together with same velocity
3. Ensure nr/qr ratio stays physical

**Implementation:**

```cpp
// WDM6 coupled sedimentation
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wdm6_sediment_rain(
    Array4<Real> const& qr_arr,
    Array4<Real> const& nr_arr,
    Array4<Real> const& den_arr,
    Array4<Real> const& delz_arr,
    Real dt,
    int i, int j, int kts, int kte
) {
    constexpr Real avtr = 841.9;
    constexpr Real bvtr = 0.8;
    
    // Top-down sedimentation
    for (int k = kte-1; k >= kts; --k) {
        if (qr_arr(i,j,k) > qcrmin && nr_arr(i,j,k) > nrmin) {
            // Calculate slope parameter (depends on both qr and nr)
            Real lamdar = wdm6_lamdar(qr_arr(i,j,k), den_arr(i,j,k), nr_arr(i,j,k));
            
            // Terminal velocity
            Real vt = avtr * std::pow(1.0 / lamdar, bvtr);
            
            // Flux out of this cell
            Real flux_qr = qr_arr(i,j,k) * vt * dt / delz_arr(i,j,k);
            Real flux_nr = nr_arr(i,j,k) * vt * dt / delz_arr(i,j,k);
            
            // Limit flux to avoid negative values
            flux_qr = amrex::min(flux_qr, qr_arr(i,j,k));
            flux_nr = amrex::min(flux_nr, nr_arr(i,j,k));
            
            // Update this cell
            qr_arr(i,j,k) -= flux_qr;
            nr_arr(i,j,k) -= flux_nr;
            
            // Add to cell below (if not at bottom)
            if (k > kts) {
                qr_arr(i,j,k-1) += flux_qr * delz_arr(i,j,k) / delz_arr(i,j,k-1);
                nr_arr(i,j,k-1) += flux_nr * delz_arr(i,j,k) / delz_arr(i,j,k-1);
            }
            
            // Accumulate surface precipitation (at k=kts)
            if (k == kts) {
                // flux_qr is the amount leaving the bottom
                // Convert to mm: (kg/kg) * (kg/m^3) * (m) = kg/m^2 = mm
                Real precip_mm = flux_qr * den_arr(i,j,k) * delz_arr(i,j,k);
                // Add to accumulation array
            }
        }
    }
}
```

**Reference:** WRF module_mp_wdm6.F lines 2700-2900 (nislfv_rain_plm6)

**Note:** WRF uses more sophisticated PLM (piecewise linear method) scheme. Above is simplified first-order upwind for initial implementation.

---

## Step 6: Integration (30 min)

### 6.1: Remove Stub Code

**Current ERF_AdvanceWDM6.cpp lines 74-256:**
- Contains stub with print statements
- Replace entirely with new implementation

### 6.2: Assembly Checklist

Verify all pieces are in place:

```cpp
// Device helper functions (lines 20-300)
[ ] wdm6_cpmcal, wdm6_xlcal, etc. (already done ✅)
[ ] wdm6_lamdar (NEW - slope with nr)
[ ] wdm6_lamdac (NEW - slope with nc)
[ ] wdm6_mean_droplet_diameter (NEW)
[ ] wdm6_ccn_activation (NEW)
[ ] wdm6_sediment_rain (NEW)

// Main Advance() function (lines 300-2500)
[ ] Physical constants setup
[ ] MFIter loop with nn, nc, nr arrays
[ ] Working FArrayBox allocations
[ ] Initialization loop (enforce minimums)
[ ] Saturation calculation (from WSM6)
[ ] CCN activation call (NEW)
[ ] Condensation/evaporation (from WSM6, minor mods)
[ ] Autoconversion (MODIFIED for nc→nr)
[ ] Accretion (from WSM6, minor mods)
[ ] Ice processes (from WSM6, no changes)
[ ] Sedimentation (MODIFIED for coupled qr/nr)
[ ] Precipitation accumulation
[ ] Final bounds enforcement
```

### 6.3: Fortran Bridge Path

Keep the `#ifdef ERF_USE_WDM6_FORT` path intact:
- Lines 95-253 currently have Fortran bridge
- Wrap new C++ code in `#else` branch
- Both paths should coexist

```cpp
#ifdef ERF_USE_WDM6_FORT
    // ... existing Fortran bridge code (lines 95-253) ...
#else
    // NEW C++ GPU kernel path goes here
    // Copy structure from WSM6
    // Add WDM6 modifications as described above
#endif
```

---

## Step 7: Testing & Validation (1-2 hours)

### 7.1: Compilation Test

```bash
cd /g/g10/wise14/compiling/clean/ERF
rm -rf build
cmake -DERF_PRECISION=DOUBLE -DERF_ENABLE_CUDA=ON -B build
cmake --build build -j8
```

**Expected errors to fix:**
- Missing variable declarations
- Array dimension mismatches
- Undefined constants (add to ERF_WDM6.H)

### 7.2: Runtime Test - No Crash

```bash
./build/erf3d Tests/WDM6_verification/inputs_bubble_wdm6
```

**Success criteria:**
- Runs without crashing
- Produces plotfiles
- No NaN values

### 7.3: Physics Validation

**Check:**
```bash
python Tests/WDM6_verification/analyze_wdm6.py plt_wdm6_bubble00500
```

**Expected results:**
- ✅ nc in range 1e1 - 1e9 m^-3
- ✅ nr in range 1e-2 - 1e6 m^-3
- ✅ nn decreases as CCN activate
- ✅ Clouds form (qc > 0)
- ✅ Rain forms (qr > 0)
- ✅ Precipitation accumulates

### 7.4: Compare with WSM6

Run same case with WSM6:
```bash
./build/erf3d Tests/WDM6_verification/inputs_bubble_wsm6
python Tests/WDM6_verification/analyze_wdm6.py plt_wdm6_bubble00500 plt_wsm6_bubble00500
```

**Expected:**
- Similar precipitation amounts (within 20%)
- Similar cloud/rain fields
- WDM6 shows more realistic droplet sizes

---

## Known Simplifications (Can Enhance Later)

### Phase 1 Implementation Uses:

1. **Simplified CCN Activation**
   - Abdul-Razzak & Ghan (2000) but simplified
   - Single aerosol mode (not multiple modes)
   - Full version in WRF lines 2040-2200

2. **First-Order Sedimentation**
   - Upwind scheme instead of PLM
   - Simpler but more diffusive
   - Full PLM in WRF lines 2700-2900

3. **Basic Collection Kernels**
   - Long's collection kernel simplified
   - Full version has more detailed size-dependent efficiency

### Enhancement Opportunities:

After basic version works, can add:
- Multi-modal aerosol activation
- PLM sedimentation scheme
- Detailed collection kernels
- Breakup parameterization
- Evaporation of rain drops

---

## File Modification Summary

### Files to Modify:

1. **ERF_AdvanceWDM6.cpp** (PRIMARY)
   - Lines 20-70: Add new device functions
   - Lines 74-256: Replace stub with full implementation
   - Estimated final size: ~2500 lines (similar to WSM6)

2. **ERF_WDM6.H** (MINOR)
   - May need to add constants (qck1, di15, etc.)
   - Already has most parameters ✅

3. **ERF_InitWDM6.cpp** (MINOR)
   - May need to add coefficient calculations
   - Check lines 30-60 for new coefficients

### Files to Reference (Don't Modify):

- `/g/g10/wise14/compiling/clean/ERF/Source/Microphysics/WSM6/ERF_AdvanceWSM6.cpp`
- `/p/lustre1/wise14/wrf/compiling/WRF/WRF_v4.7.1_clean/phys/module_mp_wdm6.F`

---

## Code Snippets for Key Sections

### Initialize Coefficients (ERF_InitWDM6.cpp)

Add to `initialize_coeffs()` around line 60:

```cpp
// WDM6-specific coefficients
m_qc0 = 4.0/3.0*M_PI*rhoh2o*std::pow(r0,3)*xncr0/rho_0;
m_qc1 = 4.0/3.0*M_PI*rhoh2o*std::pow(r0,3)*xncr1/rho_0;
m_qck1 = 0.104*CONST_GRAV*peaut/std::pow(rhoh2o, 1.0/3.0)/xmyu*std::pow(rho_0, 4.0/3.0);
m_pidnc = M_PI*rhoh2o/6.0;
m_pidnr = 4.0*M_PI*rhoh2o;

// Store in class members (add to ERF_WDM6.H):
// Real m_qc0, m_qc1, m_qck1, m_pidnc, m_pidnr;
```

### Key Constants to Add (ERF_WDM6.H)

Around line 180, add to class:

```cpp
// WDM6 double-moment parameters
Real m_qc0{0.0};      // CCN activation threshold (low)
Real m_qc1{0.0};      // CCN activation threshold (high)
Real m_qck1{0.0};     // Autoconversion coefficient
Real m_pidnc{0.0};    // pi * rho_water / 6
Real m_pidnr{0.0};    // 4 * pi * rho_water
Real m_di15{15.e-6};  // 15 micron diameter threshold
Real m_di82{82.e-6};  // 82 micron diameter (rain evap)
```

---

## Debug Checklist

When things go wrong:

### Compiler Errors

- [ ] All new variables declared in ERF_WDM6.H?
- [ ] All device functions have AMREX_GPU_HOST_DEVICE?
- [ ] No host-only functions called from device?
- [ ] Array dimensions match (i,j,k vs i,k,j)?

### Runtime Crashes

- [ ] Check for division by zero (qc=0, nc=0)
- [ ] Check array bounds (kts, kte, etc.)
- [ ] Check for NaNs with erf.check_nan = 1
- [ ] Check minimums enforced (ncmin, nrmin)

### Physics Issues

- [ ] No clouds forming: Check CCN activation, supersaturation
- [ ] Excessive clouds: Check autoconversion threshold
- [ ] No rain: Check autoconversion, sedimentation
- [ ] Unrealistic nc: Check initialization, activation rate
- [ ] nn going negative: Check activation doesn't over-deplete

---

## Success Criteria

Phase 2 is complete when:

1. ✅ Code compiles with no errors
2. ✅ Runs bubble test without crashing
3. ✅ Produces clouds (qc > 0)
4. ✅ Produces rain (qr > 0)
5. ✅ Number concentrations in physical range:
   - nc: 1e1 to 1e9 m^-3
   - nr: 1e-2 to 1e6 m^-3
   - nn: 0 to ccn0
6. ✅ Precipitation accumulates
7. ✅ Results similar to WSM6 (within factor of 2)
8. ✅ No NaN values

---

## Quick Reference: Key Line Numbers

### In WSM6 (ERF_AdvanceWSM6.cpp):

- Lines 890-960: Main loop setup
- Lines 1025-1095: FArrayBox allocations
- Lines 1146-1200: Initialization
- Lines 1600-1650: Autoconversion
- Lines 2000-2400: Sedimentation

### In WRF WDM6 (module_mp_wdm6.F):

- Lines 2085-2211: wdm6init
- Lines 334-2035: wdm62d (main physics)
- Lines 2040-2200: activ_conc (CCN activation)
- Lines 2700-2900: nislfv_rain_plm6 (sedimentation)

### In WDM6 Target (ERF_AdvanceWDM6.cpp):

- Lines 20-67: Helper functions (done ✅)
- Lines 74-256: Stub to replace
- Target size: ~2500 lines

---

## Resources

**WRF WDM6 Papers:**
- `/g/g10/wise14/compiling/clean/ERF/WDM6_paper.pdf`
- `/g/g10/wise14/compiling/clean/ERF/WDM6_description.pdf`

**Key References:**
- Abdul-Razzak & Ghan (2000): CCN activation
- Lim & Hong (2010): WDM6 scheme description

**ERF Documentation:**
- `Source/Microphysics/WDM6/README`
- `Tests/WDM6_verification/README.md`

---

## Tomorrow's Pickup

When resuming:

1. **Start here:** Step 1 - Copy WSM6 Structure
2. **First commit:** After Step 2 (number concentration tracking compiles)
3. **Second commit:** After Step 3 (CCN activation added)
4. **Third commit:** After Steps 4-5 (autoconversion + sedimentation)
5. **Test commit:** After Step 7 validation

Estimated session time: 4-6 hours with testing

Good luck! 🚀
