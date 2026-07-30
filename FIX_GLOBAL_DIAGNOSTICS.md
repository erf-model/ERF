# Fix: Global Domain Diagnostics for WDM6

## Problem Identified

Your **plotfiles showed moisture exists** (qv, qc, qp), but **console diagnostics showed zeros**. 

**Root cause:** The diagnostic code only checked the **first MFIter tile**, not the entire domain. If your supercell storm was in a different tile, the diagnostics would show zeros even though moisture and precipitation existed elsewhere.

## Changes Made

### 1. Global diagnostics BEFORE microphysics call
**Location:** `ERF_AdvanceWDM6.cpp` after line 270

**What it does:**
- Checks **entire domain** for maximum values of all moisture variables
- Runs **before** the MFIter loop, so it scans all tiles
- Uses `reduce_ops.eval(*mic_fab_vars[...])` which operates on the full MultiFab

**Output format:**
```
WDM6::Advance() call #N (dt=Xs)
  GLOBAL max mixing ratios (g/kg): qv=X, qc=X, qr=X, qi=X, qs=X, qg=X
  GLOBAL max number conc (#/kg): nc=X, nr=X, nn=X
```

### 2. Global precipitation diagnostic AFTER microphysics
**Location:** `ERF_AdvanceWDM6.cpp` before final closing brace

**What it does:**
- Checks **entire domain** for precipitation accumulation this timestep
- Runs **after** all microphysics processing complete

**Output format:**
```
  GLOBAL precip this step: rain_sum=X mm, rain_max=X mm
```

### 3. Removed misleading tile-based diagnostics
The old diagnostics that only checked one tile have been removed to avoid confusion.

## What to Look For After Rebuild

### Expected output at initialization (call #1):
```
WDM6::Advance() call #1 (dt=1s)
  GLOBAL max mixing ratios (g/kg): qv=15.5, qc=0.001, qr=0.01, qi=0.02, qs=0, qg=0
  GLOBAL max number conc (#/kg): nc=1.18e8, nr=1147, nn=1.18e9
  GLOBAL precip this step: rain_sum=0 mm, rain_max=0 mm
```

The key values to check:

1. **qv (water vapor):** Should be **10-20 g/kg** for supercell environment
   - If 0: No moisture in initialization → **fix initial conditions**
   - If >10: Good, moisture available for cloud formation

2. **qc (cloud water):** Should grow from ~0 to 1-5 g/kg as storm develops
   - If stays 0: Clouds not forming → check if saturation is reached

3. **qr (rain):** Should be 0.1-5 g/kg in mature storm
   - If stays << 0.1: Autoconversion too weak (this is your current problem!)

4. **Qp = qr + qs + qg:** Your target diagnostic
   - **Expected:** 2-20 g/kg in storm core
   - **Your current:** Order of magnitude smaller than other schemes

### Expected output during storm (call #1000):
```
WDM6::Advance() call #1000 (dt=1s)
  GLOBAL max mixing ratios (g/kg): qv=14.2, qc=2.5, qr=3.8, qi=0.5, qs=1.2, qg=2.1
  GLOBAL max number conc (#/kg): nc=2.5e8, nr=8500, nn=9.8e8
  GLOBAL precip this step: rain_sum=0.15 mm, rain_max=0.008 mm
```

## Diagnosis Based on New Output

### Case A: qv is present but qr/qs/qg are too small

**Example:**
```
qv=15 g/kg (good), qc=3 g/kg (good), but qr=0.05 g/kg (too small!)
```

**Diagnosis:** Warm rain processes (autoconversion, accretion) are too weak.

**Possible causes:**
1. **nc (droplet number) too high** → small droplets → slow autoconversion
   - Check: nc > 5e8 #/kg is very high (>500 cm^-3)
   - Fix: Reduce activation fraction in `ERF_InitWDM6.cpp`

2. **Autoconversion threshold or rate too conservative in Fortran code**
   - WDM6 uses droplet number to calculate rates
   - If nc/qc ratio is off, rates will be wrong

3. **Timestep subcycling issue**
   - Check if "Minor timesteps" is reasonable
   - Should be 1-5 for dt=1s

### Case B: qv is zero or very low

**Example:**
```
qv=0.1 g/kg (way too low!)
```

**Diagnosis:** Not enough moisture in initial conditions.

**Fix:** Check your initial sounding/profile. Supercell needs:
- Surface: qv ~ 14-18 g/kg (dewpoint ~ 20-25°C)
- 850 mb: qv ~ 10-12 g/kg
- 500 mb: qv ~ 3-5 g/kg

### Case C: Values look reasonable but Qp still too small

**Example:**
```
qv=14 g/kg, qc=2 g/kg, qr=1.5 g/kg, qs=0.8 g/kg, qg=1.2 g/kg
Qp = qr + qs + qg = 3.5 g/kg
```

**Diagnosis:** WDM6 is working, but producing less precipitation than other schemes.

**Possible reasons:**
1. **Different scheme physics:** WDM6 double-moment may distribute precipitation differently
2. **More realistic?** Some schemes over-predict precipitation
3. **Ice processes:** WDM6 ice physics may differ from scheme you're comparing to

**What to check:**
- Compare vertical profiles (not just maximum values)
- Check column-integrated values (precipitable water)
- Look at time evolution (is peak precip just delayed?)

## Comparing with Other Schemes

When you say "order of magnitude smaller," what's the comparison?

**If comparing to Morrison or P3:**
- Those are also double-moment schemes
- Should be similar
- If WDM6 is 10× smaller, something is wrong

**If comparing to single-moment schemes (Kessler, WSM6, Thompson):**
- May have different characteristic values
- Single-moment can be more aggressive
- Factor of 2-3 difference is normal; factor of 10 is not

## Next Steps

1. **Rebuild and rerun:**
   ```bash
   make clean
   make -j8
   ./ERF your_input_file
   ```

2. **Check the new GLOBAL diagnostic output** and share:
   - qv values (especially at t=0 and mature storm time)
   - qc, qr, qi, qs, qg evolution
   - Whether they match your plotfile analysis

3. **Compare maximum values to plotfile:**
   ```python
   # In your analysis script
   import yt
   ds = yt.load("plt03600")
   print(f"Max qv = {ds.r['qv'].max()*1000:.2f} g/kg")
   print(f"Max qr = {ds.r['qr'].max()*1000:.2f} g/kg")
   # Should match console "GLOBAL max mixing ratios"
   ```

4. **If values now make sense,** we can diagnose why WDM6 produces less Qp than expected:
   - Compare nc (droplet number) between schemes
   - Check if autoconversion threshold is different
   - Look at ice process rates

## What This Should Reveal

The global diagnostics will tell us:

✅ **If moisture is present** (qv > 10 g/kg) → Initial conditions are OK

✅ **If clouds are forming** (qc growing from 0 to 1-5 g/kg) → Saturation and activation working

✅ **If rain is forming** (qr > 0.1 g/kg) → Autoconversion working (even if weak)

✅ **If ice is forming** (qi, qs, qg > 0.1 g/kg above freezing level) → Ice processes working

Then we can focus on **why the rates are smaller than expected** rather than debugging whether things work at all.
