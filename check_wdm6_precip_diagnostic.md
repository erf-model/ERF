# WDM6 Low Precipitation Diagnostic Guide

## Problem
You're seeing lower Qp (= qr + qs + qg) values than expected in supercell storm core.

## What Qp represents
- **Qp = qr + qs + qg** (rain + snow + graupel mixing ratios in kg/kg)
- In a supercell updraft, you should see significant suspended precipitation
- Typical values: 1-10 g/kg in moderate updrafts, up to 20+ g/kg in strong cores

## Diagnostic Steps

### 1. Check Console Output
The code now prints diagnostic information every timestep. Look for:

```
WDM6::Advance() call #N (dt=X s)
  Number conc (#/kg): nn=[min to max], nc=[min to max], nr=[min to max]
  Max mixing ratios (g/kg): qc=X, qr=X, qi=X, qs=X, qg=X
  Precip this step: rain_sum=X mm, rain_max=X mm
```

**Expected values for a supercell:**
- `nc`: 10^8 to 10^9 #/kg (100-1000 cm^-3 activated droplets)
- `nr`: 10^3 to 10^6 #/kg (rain drops)
- `qc`: 1-5 g/kg (cloud water)
- `qr`: 0.1-10 g/kg (rain water)
- `qs`, `qg`: 0.1-5 g/kg (ice species)

### 2. Potential Issues and Fixes

#### Issue A: Low cloud droplet number (nc)
**Symptom:** `nc < 10^7 #/kg` or `nc_max < 100 cm^-3`
**Cause:** Insufficient aerosol activation
**Check:** Look at `nn` (aerosol number). Should be ~10^8-10^9 #/kg for CCN0=100e6 m^-3
**Fix location:** `ERF_InitWDM6.cpp:88-123` - activation logic

#### Issue B: Excessive sedimentation
**Symptom:** `qr, qs, qg` drop rapidly between timesteps
**Cause:** Terminal velocities too high or timestep too large
**Check WRF WDM6 terminal velocity parameters:**
- Rain: V = 841.9 × (1/λ)^0.8 m/s
- Snow: V = 11.72 × (1/λ)^0.41 m/s
- Graupel: V = 330 × (1/λ)^0.89 m/s

**Current implementation:** Simple top-down sedimentation (not PLM)
- May be too diffusive or have wrong CFL limit
- Fortran code uses full PLM scheme

#### Issue C: Autoconversion rate too slow
**Symptom:** High `qc` (cloud water) but low `qr` (rain)
**Cause:** Not converting cloud to rain fast enough
**Fix:** This is handled by Fortran `mp_wdm6_run_c()` - check if it's being called correctly

#### Issue D: Wrong physics constants
**Symptom:** All mixing ratios seem off by factor of 10 or 100
**Cause:** Unit conversion error
**Check locations:**
- Line 290-307: Physical constants passed to Fortran
- Line 276-283: WDM6 initialization constants
- Make sure `ccn0` is in correct units (m^-3, not cm^-3)

### 3. Quick Validation Checks

#### Check 1: Mass conservation
```bash
# In your analysis script, verify:
# Sum of all water species should be conserved
qtotal = qv + qc + qi + qr + qs + qg
# Should remain constant or only decrease via surface precipitation
```

#### Check 2: Physically reasonable ranges
```python
# Typical supercell values:
# qv: 5-15 g/kg (water vapor - depends on level)
# qc: 0.5-5 g/kg (cloud droplets in updraft)
# qr: 0.1-10 g/kg (rain - highest near melting level)
# qi: 0.1-2 g/kg (ice crystals)
# qs: 0.5-5 g/kg (snow/aggregates)
# qg: 1-10 g/kg (graupel/hail in strong updrafts)
# Qp = qr + qs + qg: typically 2-20 g/kg in core
```

#### Check 3: Vertical profile
In a supercell updraft, Qp should show:
- Low at surface (precipitation falls out)
- **Peak at mid-levels** (5-8 km) where conversion is active
- Moderate aloft (ice processes)

If Qp is flat or decreasing with height, sedimentation is removing it too fast.

### 4. Comparison with WRF

To validate against WRF WDM6:
1. Run same initial conditions in WRF with WDM6
2. Extract qr, qs, qg at same time
3. Compare:
   - Maximum values
   - Vertical profiles
   - Horizontal structure

### 5. Code Locations to Check

**Fortran bridge call:** `ERF_AdvanceWDM6.cpp:417-431`
- Verify array pointers are correct
- Check that dt_advance is reasonable (should be < 10 seconds for microphysics)
- Verify all constants are in correct units

**Number concentration initialization:** `ERF_InitWDM6.cpp:88-131`
- Check nn, nc, nr initialization
- Verify ccn0 / rho conversion is correct

**Sedimentation (C++ path, lines 847-948):**
- Currently simplified - may not match Fortran
- For validation, USE FORTRAN PATH (`-DERF_ENABLE_WDM6_FORT=ON`)

### 6. Recommended Next Steps

1. **Enable more diagnostics:** Add prints to show:
   - Column-integrated values
   - Maximum Qp location (i,j,k)
   - Tendencies (dqr/dt, dqs/dt, dqg/dt)

2. **Check timestep splitting:** WDM6 should sub-cycle if dt is large
   - Look for `wdm6_loops` in output
   - Should be > 1 if dt_advance > ~5-10 seconds

3. **Verify Fortran path is being used:**
   - Check for "WDM6 Fortran bridge initialized" in output
   - Confirm `ERF_USE_WDM6_FORT` is defined

4. **Compare Fortran vs C++ path:**
   - Currently C++ path has simplified ice physics
   - Fortran path has full WRF WDM6
   - Use Fortran for validation

## What to Report Back

Please share:
1. The console output showing the diagnostic prints
2. Typical values you're seeing for qr, qs, qg (g/kg)
3. What you expect to see vs what you're getting
4. Whether you're using Fortran or C++ path
5. The timestep (dt_advance) being used
6. A vertical profile plot of Qp if possible

This will help diagnose whether the issue is:
- Physics (wrong process rates)
- Numerics (sedimentation, timestep)
- Initialization (number concentrations)
- Units/constants (factor of 10 errors)
