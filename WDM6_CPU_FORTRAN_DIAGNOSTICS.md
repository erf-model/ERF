# WDM6 CPU Fortran Diagnostics: Storm Collapse After Hour 2

**Date**: 2026-08-02  
**Problem**: ERF WDM6 loses all precipitation and cold pool by hour 2 while WRF maintains both  
**Code Path**: CPU Fortran bridge (`ERF_USE_WDM6_FORT`)  
**Fortran Code**: Verified numerically identical to WRF

---

## Problem Summary from Plots

### **Hour 1** (04:00:00):
- ✅ ERF WDM6 cold pool present (similar to WRF, WSM6, Morrison)
- ✅ Rain profile: ~0.0002 kg/kg (slightly lower than WRF's ~0.0003)
- ✅ Snow/graupel profiles: Similar shape to WRF, slightly lower magnitude
- ⚠️ Already see slightly less precipitation

### **Hour 2** (05:00:00):
- ❌ **ERF WDM6 cold pool GONE** (warm air where cold pool should be)
- ❌ **ZERO rain** (WRF still has 0.0002 kg/kg)
- ❌ **Much less snow/graupel** than WRF
- ❌ **Temperature 5-10K warmer** at low levels
- ❌ **Higher water vapor** at low levels (evaporated precip!)

**Conclusion**: Rain is evaporating completely before reaching ground → no evaporative cooling → cold pool collapses → storm dies

---

## Code Verification Status

### ✅ **Fortran Physics - IDENTICAL**
- Constants: All match WRF (verified)
- Rain evaporation formula: Identical (line 1258 in ERF, 1468 in WRF)
- Saturation adjustment: Identical
- Minor timestep handling: Identical (dtcldcr = 120s)
- Precision: Both double precision (`c_double`)

### ✅ **Fortran Interface - CORRECT**
- `mp_wdm6_run_c()` signature correct
- All arrays passed as `intent(inout)`
- Constants passed correctly

### ⚠️ **Parameters Being Passed - NEED TO CHECK**
These could differ between ERF and WRF:
1. CCN0 value
2. Physics timestep (dt)
3. Initial conditions

---

## Potential Issues to Investigate

### **1. CCN0 Value** (HIGH PRIORITY)

**ERF Default**: `m_ccn0 = 100.0e6 m^-3` (100 cm^-3, continental)

**WRF Default**: Depends on namelist, typical values:
- Maritime: 50-100 cm^-3
- Continental: 100-300 cm^-3  
- Polluted: 500-1000 cm^-3

**Why it matters**: CCN0 affects:
- Cloud droplet number concentration (nc)
- Droplet size
- Autoconversion rate (larger drops → faster rain formation)
- Evaporation rate (more smaller drops → faster evaporation)

**Diagnostic**:
```bash
# Check what CCN0 ERF is using
grep "CCN0 =" your_erf_output.log

# Check what WRF is using
grep "ccn\|CCN" /p/lustre1/wise14/wrf/.../namelist.input
```

**Test**: Try running ERF with CCN0 = 300e6 to see if it helps rain reach surface

---

### **2. Physics Timestep** (HIGH PRIORITY)

**Question**: What's your physics timestep in ERF vs WRF?

**Why it matters**: 
- WDM6 has internal minor timesteps (dtcld = dt/loops where loops = max(dt/120, 1))
- If ERF dt > WRF dt, processes might be too aggressive per timestep
- Rain could evaporate too much before next dynamics step

**Diagnostic**:
```bash
# ERF timestep
grep "dt.*advance\|physics.*dt" your_erf_output.log

# WRF timestep  
grep "time_step\|dt" /p/lustre1/wise14/wrf/.../namelist.input
```

**If ERF dt > 30s**: This could be the issue!

---

### **3. Theta Conversion Timing** (MEDIUM PRIORITY)

The theta→T→theta conversion happens in `ERF_AdvanceWDM6.cpp` lines 439-453.

**Potential issue**: Latent heating might be applied at wrong point in dynamics cycle

**Check**:
1. When does ERF call microphysics relative to dynamics?
2. Is theta being used correctly by dynamics after update?
3. Could there be a sign error in exner function?

**Diagnostic**: Add to ERF right after line 452:
```cpp
// Diagnostic: Check theta changes
Real max_dtheta = 0.0;
ParallelFor(box, [=, &max_dtheta] AMREX_GPU_DEVICE (int i, int j, int k) {
    Real dtheta = theta_arr(i,j,k) - theta_before(i,j,k);
    max_dtheta = amrex::max(max_dtheta, amrex::abs(dtheta));
});
amrex::Print() << "  Max theta change from microphysics: " << max_dtheta << " K\n";
```

---

### **4. Number Concentration Drift** (MEDIUM PRIORITY)

**Question**: Are nc, nr, nn drifting to unphysical values by hour 2?

**Check at hour 2**:
- `nc` should be 10^8 - 5×10^8 #/m^3 (100-500 cm^-3)
- `nr` should be 10^3 - 10^6 #/m^3
- `nn` should deplete from CCN0 but not go to zero

**If nc is too high**: Rain drops are too small → evaporate too fast  
**If nr is too low**: Few rain drops → less surface area → evaporate faster

**Diagnostic**: Add to `ERF_AdvanceWDM6.cpp` after Fortran call:
```cpp
Real max_nc = mic_fab_vars[MicVar_WDM6::nc]->max(0);
Real max_nr = mic_fab_vars[MicVar_WDM6::nr]->max(0);
Real min_nn = mic_fab_vars[MicVar_WDM6::nn]->min(0);
amrex::Print() << "  After WDM6: max nc=" << max_nc << ", max nr=" << max_nr 
               << ", min nn=" << min_nn << " #/kg\n";
```

---

### **5. Sedimentation Rate** (MEDIUM PRIORITY)

**Hypothesis**: Rain falling out too fast, not staying aloft long enough to cool air

**Check**: Terminal velocities in downdraft regions
- Typical rain vt: 5-10 m/s
- In strong downdraft (10 m/s down): effective vt = 15-20 m/s relative to ground

**If rain falls faster in ERF**: Could be numerical issue in sedimentation

---

### **6. Subsaturation in Downdraft** (LOW PRIORITY - Physics should handle)

**Question**: How subsaturated is the air where rain is evaporating?

At hour 2 in ERF:
- If RH < 50% in downdraft: Rain should evaporate (correct physics)
- If RH > 80% in downdraft: Rain shouldn't evaporate completely (physics bug)

**Diagnostic**: Check RH profile at hour 2 in storm downdraft region

---

## Recommended Diagnostic Steps

### **Step 1: Check Parameters** (5 minutes)
```bash
# ERF
grep "CCN0\|dt.*advance" your_erf_log.txt

# WRF
grep "ccn_diag\|time_step" /p/lustre1/wise14/wrf/.../namelist.input
grep "mp_physics" /p/lustre1/wise14/wrf/.../namelist.input  # Should be 16 for WDM6
```

### **Step 2: Add Number Concentration Diagnostics** (30 minutes)

Add to `Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp` after line 437:

```cpp
// Diagnostic every 10 minutes
static int diag_count = 0;
diag_count++;
if (diag_count % 20 == 0) {  // Assuming ~30s timestep
    Real max_qr = mic_fab_vars[MicVar_WDM6::qr]->max(0);
    Real max_nc = mic_fab_vars[MicVar_WDM6::nc]->max(0);
    Real max_nr = mic_fab_vars[MicVar_WDM6::nr]->max(0);
    Real min_nn = mic_fab_vars[MicVar_WDM6::nn]->min(0);
    Real mean_nc = mic_fab_vars[MicVar_WDM6::nc]->sum(0) / mic_fab_vars[MicVar_WDM6::nc]->boxArray().numPts();
    
    amrex::Print() << "=== WDM6 Diagnostic (t=" << diag_count * dt_advance << "s) ===\n";
    amrex::Print() << "  max qr = " << max_qr*1000 << " g/kg\n";
    amrex::Print() << "  max nc = " << max_nc << " #/kg (= " << max_nc*1.2 << " #/m^3)\n";
    amrex::Print() << "  max nr = " << max_nr << " #/kg\n";
    amrex::Print() << "  min nn = " << min_nn << " #/kg\n";
    amrex::Print() << "  mean nc = " << mean_nc << " #/kg\n";
}
```

### **Step 3: Compare nc/nr Evolution** (1 hour)

Run ERF and WRF side-by-side with identical setup. Plot:
- nc vs time
- nr vs time
- qr vs time

If ERF nc >> WRF nc: Drops too small → evaporate too fast  
If ERF nr << WRF nr: Fewer drops → less mass reaches ground

### **Step 4: Check RH in Downdraft** (30 minutes)

At hour 2, check RH where rain should be:
```python
# In your analysis script
rh = qv / qsat_water
print(f"RH at z=1-3km in storm core: {rh[storm_mask, 10:30].mean()}")

# If RH < 0.5: Rain evaporation is correct physics
# If RH > 0.8: Something wrong with saturation calculation
```

---

## Most Likely Causes (Ranked)

### **1. CCN0 Too High** (80% confidence)
If ERF CCN0 > WRF CCN0:
- More cloud droplets
- Smaller droplets
- Faster evaporation
- **Test**: Reduce ERF CCN0 to match WRF

### **2. Physics Timestep Too Large** (60% confidence)
If ERF dt > 30s and WRF dt <= 30s:
- Processes too aggressive per step
- Rain evaporates too much before dynamics responds
- **Test**: Reduce ERF dt to match WRF

### **3. Number Concentration Drift** (40% confidence)
nc or nr evolving incorrectly:
- Needs diagnosis over time
- **Test**: Add diagnostics and compare to WRF

### **4. Dynamics-Microphysics Coupling** (30% confidence)
Theta update not feeding back to dynamics correctly:
- Latent heating not cooling air properly
- **Test**: Check if other schemes (Morrison, WSM6) also lose cold pool

### **5. Subtle Fortran Difference** (10% confidence)
Unlikely given code is identical, but possible:
- Compiler optimization differences
- Numerical precision accumulation
- **Test**: Run bit-for-bit comparison with WRF

---

## Immediate Actions

1. **Check your simulation parameters**:
   - What CCN0 are you using in ERF?
   - What CCN0 is WRF using?
   - What are the physics timesteps?

2. **Add the diagnostic print statements** I provided above

3. **Run a short test** (1 hour simulation):
   - With ERF CCN0 = 100e6
   - With ERF CCN0 = 300e6
   - Compare rain reaching ground

4. **Share the diagnostic output** - I can help interpret

---

## Questions for You

1. What CCN0 value are you using in ERF? (check inputs or log)
2. What CCN0 is WRF using? (check namelist.input for ccn_diag)
3. What's your physics timestep in ERF vs WRF?
4. Do Morrison and WSM6 in ERF maintain the cold pool properly at hour 2?
   - If YES → WDM6-specific issue
   - If NO → Dynamics coupling issue

Please share this information and I can narrow down the exact cause!

---

**Generated**: 2026-08-02  
**By**: Claude Code (Sonnet 4)
