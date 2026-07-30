# Most Likely Issues Causing Low Qp in WDM6 Supercell

## Priority 1: Number Concentration Issues

### Issue: Cloud droplet number (nc) initialization
**Location:** `ERF_InitWDM6.cpp:88-131`

The current code initializes `nc` from aerosols, but there may be issues:

```cpp
// Line 98-122: Current initialization
if (nc(i,j,k) < Real(1.e1)) {
    if (qc(i,j,k) > Real(1.e-9)) {
        // Cloud exists - activate 15% of aerosols
        Real activation_fraction = Real(0.15);
        Real nc_init = nn(i,j,k) * activation_fraction;
        ...
```

**Problem:** If `nn` (aerosol) is initialized but aerosols aren't being activated properly, nc could be too low or too high.

**Quick Test:** Add diagnostic to check actual nc values:
- Too low nc → large droplets → fast autoconversion → precipitation falls out immediately
- Too high nc → small droplets → slow autoconversion → no rain forms

**Expected values for supercell:**
- Continental: nc ~ 200-500 cm^-3 = 2e8 - 5e8 #/m^3
- At ρ = 1.0 kg/m^3: nc ~ 2e8 - 5e8 #/kg

**Fix to try:**
```cpp
// More aggressive activation for cloud
if (qc(i,j,k) > Real(1.e-9)) {
    // Use higher activation fraction
    Real activation_fraction = Real(0.3);  // 30% instead of 15%
    ...
}
```

## Priority 2: Timestep Issues

### Issue: Microphysics timestep too large
**Location:** How `dt_advance` is set in the main code

WDM6 is sensitive to timestep. Check:
1. What is `dt_advance` in your simulation?
2. Is WDM6 sub-cycling properly?

**Expected behavior:**
- `dt_advance` should be 1-10 seconds for microphysics
- If larger, WDM6 should sub-cycle (split into smaller steps)
- Look for "Minor timesteps = N" in output

**Problem symptoms:**
- If dt is too large (> 30 sec), processes can overshoot
- Sedimentation CFL can be violated → precipitation falls through multiple cells per step
- Autoconversion/accretion can be unstable

**Quick check:** Look at your console output for:
```
WDM6::Advance() call #N (dt=X s)
```
If X > 20 seconds, this could be the issue.

## Priority 3: Sedimentation Removing Precipitation Too Fast

### Issue: Terminal velocity or sedimentation timestep
**Location:** Fortran code handles this (you're using Fortran path)

The WRF WDM6 uses PLM (piecewise linear) sedimentation which should be accurate, but:

**Potential problems:**
1. **CFL condition:** If (V_terminal * dt / dz) > 1, precipitation can fall through multiple cells
2. **Terminal velocity formulas:** Should match WRF exactly

**Expected rain terminal velocity:**
- Small drops (D ~ 0.1 mm): V ~ 0.5 m/s
- Medium drops (D ~ 1 mm): V ~ 4 m/s  
- Large drops (D ~ 5 mm): V ~ 9 m/s

**Quick check:**
- Calculate: V_max * dt / dz
- If this ratio > 0.5, you may have CFL issues
- Typical: V_max ~ 10 m/s, dz ~ 500 m, dt ~ 10 s → ratio = 0.2 (OK)
- If dz ~ 100 m, dt ~ 10 s → ratio = 1.0 (BAD!)

## Priority 4: Missing Physical Constants or Unit Errors

### Issue: CCN0 or other constants in wrong units
**Location:** `ERF_AdvanceWDM6.cpp:276-307`

Check that constants are correct:

```cpp
constexpr double den0 = 1.28;                      // kg/m^3 ✓
constexpr double denr = static_cast<double>(rhoh2o);  // Should be 1000 kg/m^3
const double ccn0 = static_cast<double>(m_ccn0);   // Should be m^-3, not cm^-3
```

**Common error:** 
- If ccn0 = 100 (thinking cm^-3) instead of 100e6 (m^-3), you get 1 million times too few aerosols!
- Check: m_ccn0 should be ~ 100e6 to 1000e6 for continental air

**Quick check in your input file:**
```
wdm6.ccn0 = 100.0e6  # Correct: 100/cm^3 = 100e6/m^3
# NOT
wdm6.ccn0 = 100.0    # Wrong! This would be 100/m^3
```

## Priority 5: Warm Rain Processes Too Weak

### Issue: Autoconversion or accretion rates
**Location:** Fortran code `ERF_module_mp_wdm6.F90`

WDM6 warm rain processes:
1. **Autoconversion:** qc → qr (cloud droplets collide to form rain)
2. **Accretion:** qc + qr → qr (rain collects cloud droplets)

These are controlled by nc (cloud droplet number).

**Relationship:**
- High nc → small droplets → slow autoconversion → weak warm rain
- Low nc → large droplets → fast autoconversion → strong warm rain

**Typical rates for supercell:**
- Autoconversion threshold: qc > 0.5 g/kg
- Autoconversion rate: ~ 0.1-1 g/kg per minute
- Accretion rate: ~ 1-10 g/kg per minute (faster)

**Diagnostic:** Check if qc builds up without forming qr
- If qc > 5 g/kg but qr < 0.1 g/kg → autoconversion is too slow
- This points to nc being too high

## Diagnostic Priority Order

Run these checks in order:

### 1. Check console output (already implemented)
Look at the diagnostic prints for typical values:
```bash
grep "Number conc\|Max mixing ratios\|Precip this step" your_output.log
```

Expected in supercell core:
- `nc`: 1e8 - 1e9 #/kg
- `qc`: 1-5 g/kg
- `qr`: 0.5-5 g/kg
- `qs + qg`: 1-10 g/kg

### 2. Check timestep
```bash
grep "dt=" your_output.log | head -1
```
Should be < 20 seconds for microphysics stability.

### 3. Check grid spacing
```bash
grep "dx\|dy\|dz" your_input_file
```
- If dz < 200 m, may need smaller dt
- If dz > 1000 m, may not resolve updraft properly

### 4. Check CCN0 setting
```bash
grep "ccn0\|CCN0" your_input_file your_output.log
```
Should see: CCN0 = 1.0e8 or 1.0e9 m^-3

### 5. Check vertical profile
Plot Qp vs height in storm core:
- Should peak at mid-levels (5-8 km)
- Should NOT be monotonically decreasing with height
- If it is decreasing, sedimentation is too strong

## Recommended Fixes to Try

### Fix 1: Increase activation efficiency (if nc is low)
Edit `ERF_InitWDM6.cpp` line 102:
```cpp
Real activation_fraction = Real(0.3);  // was 0.15
```

### Fix 2: Check CCN0 value
Edit your input file to ensure:
```
wdm6.ccn0 = 100.0e6  # 100 per cm^3 in SI units
```

### Fix 3: Reduce timestep (if dt > 20 s)
In your input, try halving the microphysics timestep.

### Fix 4: Add more diagnostics
Add this to `ERF_AdvanceWDM6.cpp` after line 473 to see process rates:

```cpp
// NEW: Check if autoconversion is working
Real qc_sum = 0.0, qr_sum = 0.0;
ReduceOps<ReduceOpSum, ReduceOpSum> sum_ops;
ReduceData<Real, Real> sum_data(sum_ops);
sum_ops.eval(box, sum_data,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        return makeTuple(qc_arr(i,j,k), qr_arr(i,j,k));
    });
auto sums = sum_data.value();
amrex::Print() << "  Domain totals: qc_sum=" << get<0>(sums) 
              << " qr_sum=" << get<1>(sums) << "\n";
// If qc_sum is large but qr_sum is small, autoconversion is weak
```

## What to Share for Further Diagnosis

Please provide:
1. **Console output excerpt** showing diagnostic values
2. **Your input file** (or relevant sections with ccn0, dt, dx, dy, dz)
3. **Vertical profile** of Qp in storm core if you can plot it
4. **Time evolution** of max(Qp) in domain - does it grow then fall, or stay flat?
5. **Comparison values** - what Qp magnitude do you expect vs what you see?

With this information, I can pinpoint the exact issue.
