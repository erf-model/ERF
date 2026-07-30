# Compilation Fix Applied

## Problem
The global diagnostic code was using the wrong AMReX reduce API signature for MultiFab operations, causing compilation errors.

## Solution
Replaced complex reduce operations with simple, built-in MultiFab functions:
- `MultiFab::max(component)` - returns maximum value across all tiles/processes
- `MultiFab::sum(component)` - returns sum across all tiles/processes

## Changes Made

### ERF_AdvanceWDM6.cpp
**Before (complex, didn't compile):**
```cpp
ReduceOps<...> reduce_ops;
ReduceData<...> reduce_data(reduce_ops);
reduce_ops.eval(*mic_fab_vars[MicVar_WDM6::qv], reduce_data, [=] {...});
```

**After (simple, works):**
```cpp
Real max_qv = mic_fab_vars[MicVar_WDM6::qv]->max(0);
Real max_qc = mic_fab_vars[MicVar_WDM6::qc]->max(0);
// etc for all variables
```

### ERF_InitWDM6.cpp
**Changed:** Loop over MFIter tiles, use local reduce per tile, then reduce across MPI ranks

## What the Diagnostics Do

### 1. Copy_State_to_Micro (first call only)
```
Copy_State_to_Micro diagnostic:
  Max qv in state = X g/kg
  Max qc in state = X g/kg
```
Shows if moisture exists in the ERF state **before** copying to microphysics.

### 2. WDM6::Advance (every call)
```
WDM6::Advance() call #N (dt=Xs)
  GLOBAL max mixing ratios (g/kg): qv=X, qc=X, qr=X, qi=X, qs=X, qg=X
  GLOBAL max number conc (#/kg): nc=X, nr=X, nn=X
```
Shows maximum values across **entire domain** (all tiles, all MPI ranks).

### 3. Precipitation (every call)
```
  GLOBAL precip this step: rain_sum=X mm, rain_max=X mm
```
Shows total and maximum precipitation this timestep across entire domain.

## Expected Compilation
```bash
make -j8
```
Should now compile successfully without errors.

## Expected Output

### At t=0 (call #1):
```
Copy_State_to_Micro diagnostic:
  Max qv in state = 15.5 g/kg
  Max qc in state = 0.0 g/kg

WDM6::Advance() call #1 (dt=1s)
  GLOBAL max mixing ratios (g/kg): qv=15.5, qc=0.001, qr=0.01, qi=0.02, qs=0, qg=0
  GLOBAL max number conc (#/kg): nc=1.18e8, nr=1147, nn=1.18e9
  GLOBAL precip this step: rain_sum=0 mm, rain_max=0 mm
```

### During storm development (call #1000):
```
WDM6::Advance() call #1000 (dt=1s)
  GLOBAL max mixing ratios (g/kg): qv=14.2, qc=2.5, qr=3.8, qi=0.5, qs=1.2, qg=2.1
  GLOBAL max number conc (#/kg): nc=2.5e8, nr=8500, nn=9.8e8
  GLOBAL precip this step: rain_sum=0.15 mm, rain_max=0.008 mm
```

## What to Check

1. **Does qv match your plotfile?**
   - Console should show qv ~ 10-20 g/kg if moisture is present
   - Compare to plotfile maximum

2. **Does qr match your plotfile?**
   - This is the key - you said plotfiles show "very small" qr
   - Console should now confirm this globally

3. **Is qc (cloud water) present?**
   - If qv > 10 but qc = 0, clouds aren't forming
   - If qc > 0 but qr << qc, autoconversion is weak

4. **How does Qp = qr + qs + qg evolve?**
   - Track it over time
   - Compare to other microphysics schemes

## What This Tells Us

The old diagnostics said "all zeros" but plotfiles showed moisture. Now:

- **If diagnostics still show zeros:** Bug in the diagnostic code (unlikely after this fix)
- **If diagnostics match plotfiles:** Confirms WDM6 is producing less precipitation than expected
  - Then we can investigate **why**: nc too high? Autoconversion too weak? Different physics?

The global diagnostics will finally give us accurate values to diagnose the low precipitation issue!
