# WDM6 Diagnostic Comparison Guide

**Date**: 2026-07-31  
**Purpose**: Compare ERF and WRF WDM6 number concentration evolution

---

## What Was Added

### **ERF Diagnostics** (ERF_UpdateWDM6.cpp)

Added print statements after `Copy_Micro_to_State` to track max nn, nc, nr for first 10 timesteps:

```cpp
if (update_call_count <= 10) {
    Real max_nn = mic_fab_vars[MicVar_WDM6::nn]->max(0);
    Real max_nc = mic_fab_vars[MicVar_WDM6::nc]->max(0);
    Real max_nr = mic_fab_vars[MicVar_WDM6::nr]->max(0);
    
    amrex::Print() << "ERF WDM6 timestep " << update_call_count
                  << " max nn,nc,nr (#/kg): "
                  << max_nn << " " << max_nc << " " << max_nr << "\n";
}
```

**Output format**:
```
ERF WDM6 timestep 1 max nn,nc,nr (#/kg): 1.2345e+02 5.6789e+01 3.4567e-01
ERF WDM6 timestep 2 max nn,nc,nr (#/kg): 1.2340e+02 5.8900e+01 3.5600e-01
...
```

### **WRF Diagnostics** (module_mp_wdm6.F)

Added two print statements:

1. **After nn initialization** (line ~228):
   ```fortran
   if (itimestep .eq. 1) then
     ! ... initialize nn = ccn0 ...
     write(*,*) 'WRF WDM6: itimestep=1 initialized nn to ccn0=', ccn0
   endif
   ```

2. **After physics call** (line ~283):
   ```fortran
   if (itimestep .le. 10 .and. j .eq. jts) then
     max_nn_local = maxval(nn(its:ite,kts:kte,j))
     max_nc_local = maxval(nc(its:ite,kts:kte,j))
     max_nr_local = maxval(nr(its:ite,kts:kte,j))
     write(*,'(A,I3,A,3E12.4)') 'WRF WDM6 timestep ', itimestep, &
                                ' max nn,nc,nr (#/kg): ', &
                                max_nn_local, max_nc_local, max_nr_local
   endif
   ```

**Output format**:
```
WRF WDM6: itimestep=1 initialized nn to ccn0=   1.000000E+08
WRF WDM6 timestep   1 max nn,nc,nr (#/kg):   1.2345E+02  5.6789E+01  3.4567E-01
WRF WDM6 timestep   2 max nn,nc,nr (#/kg):   1.2340E+02  5.8900E+01  3.5600E-01
...
```

---

## Files Modified

### **ERF**:
- `Source/Microphysics/WDM6/ERF_UpdateWDM6.cpp` (added diagnostics)
- `Source/Microphysics/WDM6/ERF_InitWDM6.cpp` (fixes from earlier)
- `Source/Microphysics/WDM6/ERF_WDM6.H` (m_nn_initialized flag)

### **WRF**:
- `/p/lustre1/wise14/wrf/compiling/WRF/WRF_v4.7.1_WDM6/phys/module_mp_wdm6.F`

---

## How to Build and Run

### **Build ERF**:
```bash
cd /g/g10/wise14/compiling/clean/ERF
cmake --build build -j8
```

### **Build WRF**:
```bash
cd /p/lustre1/wise14/wrf/compiling/WRF/WRF_v4.7.1_WDM6
./clean -a
./compile em_real >& compile.log
```

### **Run ERF**:
```bash
cd /path/to/your/test/case
./erf3d inputs > erf_output.log 2>&1
```

### **Run WRF**:
```bash
cd /path/to/your/wrf/run
mpirun -np 4 ./wrf.exe > wrf_output.log 2>&1
```

---

## What to Compare

### **1. Initial NN Value (Timestep 1)**

**ERF**: Should see nn initialized to `ccn0/ρ`
- If ccn0 = 100e6 #/m³ and ρ ≈ 1.2 kg/m³
- Expected nn ≈ 83e6 #/kg ≈ 8.3e7 #/kg

**WRF**: Should see nn initialized to `ccn0`
- **CRITICAL**: Check what units WRF uses!
- If WRF output shows `ccn0 = 1.0E+08`, check if that's #/m³ or #/kg

**⚠️ UNITS ISSUE**:
The WRF printout will tell us definitively if ccn0 is in #/m³ or #/kg:
- If WRF prints `ccn0 = 1.0E+08` and nn ≈ 1.0E+08, then ccn0 is already in #/kg
- If WRF prints `ccn0 = 1.0E+08` and nn ≈ 8.3E+07, then ccn0 is in #/m³

### **2. Number Concentration Evolution (Timesteps 1-10)**

Look for these patterns:

#### **NN (Aerosol Number)**:
- ✅ **Should DECREASE** over time (aerosols activate to form cloud droplets)
- ✅ **ERF and WRF should track together**
- ❌ If ERF nn stays constant → re-initialization bug (should be fixed now!)

#### **NC (Cloud Droplet Number)**:
- ✅ **Should INCREASE** over time (CCN activation)
- ✅ **Should start LOW** (~10 #/kg = ncmin)
- ✅ **Build gradually** as supersaturation activates aerosols
- ❌ If ERF nc jumps to ~5e7 immediately → forced initialization bug (should be fixed now!)

#### **NR (Rain Drop Number)**:
- ✅ **Should be SMALL or ZERO initially** (no rain yet)
- ✅ **Increases only when rain forms** (autoconversion)
- ✅ **Much smaller than nc** (typically nr << nc)

### **3. Physical Ranges (Sanity Check)**

All values in **#/kg**:

| Variable | Typical Range | Physical Meaning |
|----------|---------------|------------------|
| **nn** | 1e7 - 1e8 | Aerosols: 10-100 million/kg |
| **nc** | 1e1 - 1e6 | Cloud droplets: 10 - 1 million/kg |
| **nr** | 1e-2 - 1e4 | Rain: 0.01 - 10,000/kg |

**Converting to #/m³** (multiply by ρ ≈ 1.2 kg/m³):
- nn: 1.2e7 - 1.2e8 #/m³ (12 - 120 million/m³)
- nc: 12 - 1.2e6 #/m³ (0.012 - 1200 cm⁻³)
- nr: 0.012 - 12,000 #/m³

---

## Example Comparison

### **Good Agreement** ✅:

```
ERF WDM6 timestep 1 max nn,nc,nr (#/kg): 8.3000e+07 1.0000e+01 0.0000e+00
WRF WDM6 timestep 1 max nn,nc,nr (#/kg): 8.3000E+07 1.0000E+01 0.0000E+00

ERF WDM6 timestep 2 max nn,nc,nr (#/kg): 8.2950e+07 2.5000e+02 0.0000e+00
WRF WDM6 timestep 2 max nn,nc,nr (#/kg): 8.2948E+07 2.5012E+02 0.0000E+00

ERF WDM6 timestep 5 max nn,nc,nr (#/kg): 8.2700e+07 5.8000e+04 1.2000e-01
WRF WDM6 timestep 5 max nn,nc,nr (#/kg): 8.2698E+07 5.8023E+04 1.2010E-01
```

**Interpretation**:
- ✅ nn decreasing (aerosols depleting)
- ✅ nc increasing dramatically (clouds forming)
- ✅ nr appearing (rain starting)
- ✅ Values match within ~0.1%

### **Bad Agreement** ❌ (Before Fixes):

```
ERF WDM6 timestep 1 max nn,nc,nr (#/kg): 8.3000e+07 5.0000e+07 0.0000e+00  ← nc too high!
WRF WDM6 timestep 1 max nn,nc,nr (#/kg): 8.3000E+07 1.0000E+01 0.0000E+00

ERF WDM6 timestep 2 max nn,nc,nr (#/kg): 8.3000e+07 5.2000e+07 0.0000e+00  ← nn not depleting!
WRF WDM6 timestep 2 max nn,nc,nr (#/kg): 8.2948E+07 2.5012E+02 0.0000E+00
```

**Problems**:
- ❌ ERF nc = 5e7 on first step (forced initialization)
- ❌ ERF nn not changing (replenished every step)
- ❌ ERF nc not growing naturally (already at max)

---

## Troubleshooting

### **If ERF and WRF nn differ by factor of ρ**:

**Symptom**: 
```
ERF: nn = 8.3e7 #/kg
WRF: nn = 1.0e8 #/kg  (ratio ≈ 1.2 = ρ)
```

**Problem**: Units mismatch in ccn0
- ERF assumes ccn0 in #/m³, divides by ρ
- WRF passes ccn0 already in #/kg

**Fix**: In `ERF_InitWDM6.cpp`, change:
```cpp
// OLD:
nn(i,j,k) = ccn0_local / rho(i,j,k);  // Assumes ccn0 in #/m³

// NEW:
nn(i,j,k) = ccn0_local;  // ccn0 already in #/kg
```

### **If nc jumps immediately to 50 cm⁻³**:

**Problem**: Forced initialization not removed

**Check**: Make sure fix #2 was applied correctly

### **If nn stays constant in ERF**:

**Problem**: Still re-initializing every timestep

**Check**: Make sure `m_nn_initialized` flag is working

### **If values differ by orders of magnitude**:

**Check**:
1. Same initial conditions (qv, T, P)?
2. Same ccn0 value in both codes?
3. Same domain size and resolution?
4. Same timestep?
5. Is xland the same (land vs water)?

---

## Key Questions to Answer

After running both codes, you can answer:

1. **What is ccn0 in WRF?** (Check itimestep=1 printout)
   - If 1.0e8 #/m³ → ERF is correct to divide by ρ
   - If 1.0e8 #/kg → ERF should NOT divide by ρ

2. **Does nn deplete in both codes?** (Compare timesteps 1-10)
   - If yes → Aerosol-cloud interaction working
   - If ERF constant but WRF decreases → ERF re-init bug

3. **Does nc build naturally in ERF?** (Start low, grow gradually)
   - If yes → CCN activation working correctly
   - If starts at 50 cm⁻³ → Forced initialization bug

4. **Do values agree within ~1%?** (Accounting for numerics)
   - If yes → ERF = WRF, coupling is perfect!
   - If no → Investigate specific differences

---

## Success Criteria

**After fixes, you should see**:

✅ ERF and WRF nn match within 1%  
✅ ERF and WRF nc match within 5% (more variable)  
✅ ERF and WRF nr match within 10% (even more variable)  
✅ nn decreases over time in both codes  
✅ nc starts low and increases in both codes  
✅ nr appears when rain forms  

**If all criteria met → ERF WDM6 implementation is correct!**

---

## Next Steps After Verification

1. **If values match**: 
   - Remove diagnostic prints (or reduce to timestep 1 only)
   - Commit the fixes
   - Run production cases

2. **If values differ**: 
   - Document specific differences
   - Check for remaining bugs
   - May need to investigate specific physics processes

3. **If units mismatch found**:
   - Fix ERF's nn initialization
   - Re-test
   - Document the correction

---

**Created**: 2026-07-31  
**By**: Claude Code (Sonnet 4)  
**Purpose**: Guide for comparing ERF and WRF WDM6 behavior
