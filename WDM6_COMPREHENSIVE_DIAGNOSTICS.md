# WDM6 Comprehensive Hydrometeor Diagnostics

**Date**: 2026-07-31  
**Purpose**: Compare all hydrometeor statistics between ERF and WRF WDM6

---

## What Was Added

Both ERF and WRF now output **comprehensive statistics** for all hydrometeors for the first 10 timesteps:

### Statistics Computed:
- **min**: Minimum value across domain
- **max**: Maximum value across domain  
- **mean**: Domain average

### Variables Tracked:

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| **qv** | Water vapor mixing ratio | g/kg | Always present |
| **qc** | Cloud water mixing ratio | g/kg | Liquid cloud |
| **qi** | Ice mixing ratio | g/kg | Ice crystals |
| **qr** | Rain mixing ratio | g/kg | Precipitating liquid |
| **qs** | Snow mixing ratio | g/kg | Precipitating ice |
| **qg** | Graupel mixing ratio | g/kg | Rimed ice/hail |
| **Qp** | Total precip (qr+qs+qg) | g/kg | **ERF-specific composite** |
| **nn** | Aerosol number conc | #/kg | WDM6-specific |
| **nc** | Cloud droplet number | #/kg | Double-moment |
| **nr** | Rain drop number | #/kg | Double-moment |

---

## Output Format

### **ERF Output** (printed after Copy_Micro_to_State):

```
=== ERF WDM6 Statistics: Timestep 1 ===
qv (kg/kg):  min=5.0000e-03  max=1.2000e-02  mean=8.5000e-03 g/kg
qc (kg/kg):  min=0.0000e+00  max=5.0000e-04  mean=1.2000e-05 g/kg
qi (kg/kg):  min=0.0000e+00  max=0.0000e+00  mean=0.0000e+00 g/kg
qr (kg/kg):  min=0.0000e+00  max=0.0000e+00  mean=0.0000e+00 g/kg
qs (kg/kg):  min=0.0000e+00  max=0.0000e+00  mean=0.0000e+00 g/kg
qg (kg/kg):  min=0.0000e+00  max=0.0000e+00  mean=0.0000e+00 g/kg
Qp (kg/kg):  min=0.0000e+00  max=0.0000e+00  mean=0.0000e+00 g/kg  (qr+qs+qg)
nn (#/kg):   min=8.2500e+07  max=8.3500e+07  mean=8.3000e+07
nc (#/kg):   min=1.0000e+01  max=1.0000e+01  mean=1.0000e+01
nr (#/kg):   min=1.0000e-02  max=1.0000e-02  mean=1.0000e-02
========================================
```

### **WRF Output** (printed after wdm62D call):

```
=== WRF WDM6 Statistics: Timestep   1 ===
qv (kg/kg):  min=  5.0000E-03  max=  1.2000E-02  mean=  8.5000E-03 g/kg
qc (kg/kg):  min=  0.0000E+00  max=  5.0000E-04  mean=  1.2000E-05 g/kg
qi (kg/kg):  min=  0.0000E+00  max=  0.0000E+00  mean=  0.0000E+00 g/kg
qr (kg/kg):  min=  0.0000E+00  max=  0.0000E+00  mean=  0.0000E+00 g/kg
qs (kg/kg):  min=  0.0000E+00  max=  0.0000E+00  mean=  0.0000E+00 g/kg
qg (kg/kg):  min=  0.0000E+00  max=  0.0000E+00  mean=  0.0000E+00 g/kg
Qp (kg/kg):  min=  0.0000E+00  max=  0.0000E+00  mean=  0.0000E+00 g/kg  (qr+qs+qg)
nn (#/kg):   min=  8.2500E+07  max=  8.3500E+07  mean=  8.3000E+07
nc (#/kg):   min=  1.0000E+01  max=  1.0000E+01  mean=  1.0000E+01
nr (#/kg):   min=  1.0000E-02  max=  1.0000E-02  mean=  1.0000E-02
========================================
```

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
tail compile.log  # Check for SUCCESS
```

### **Run and Capture Output**:
```bash
# ERF
cd /your/test/case
./erf3d inputs | tee erf_wdm6_diag.log

# WRF  
cd /your/wrf/run
mpirun -np 4 ./wrf.exe | tee wrf_wdm6_diag.log
```

### **Extract Statistics**:
```bash
# ERF
grep -A 12 "ERF WDM6 Statistics" erf_wdm6_diag.log > erf_stats.txt

# WRF
grep -A 12 "WRF WDM6 Statistics" wrf_wdm6_diag.log > wrf_stats.txt

# Side-by-side comparison
paste erf_stats.txt wrf_stats.txt | less
```

---

## What to Look For

### **1. Initial Conditions (Timestep 1)**

#### **Water Vapor (qv)**:
- ✅ Should be **same in both codes** (initial condition)
- ✅ Typical range: 1-15 g/kg depending on temperature/altitude
- ✅ Should NOT be zero anywhere (unless very cold/dry)

#### **Cloud Water (qc)**:
- ✅ If warm start: may have some qc initially
- ✅ If cold start: should be zero, build gradually
- ⚠️ **Check max qc**: Should NOT jump to large values immediately

#### **Number Concentrations**:
- ✅ **nn**: Should be ~ccn0/ρ ≈ 8.3e7 #/kg (if ccn0=100e6, ρ=1.2)
- ✅ **nc**: Should be ~10 #/kg (ncmin) if cold start
- ✅ **nr**: Should be ~0.01 #/kg (nrmin)

### **2. Evolution (Timesteps 1-10)**

#### **Aerosol Depletion (nn)**:
```
Expected pattern:
Step 1:  nn ≈ 8.3e7  (initialized)
Step 2:  nn ≈ 8.29e7 (slightly decreasing)
Step 5:  nn ≈ 8.2e7  (continuing to decrease)
Step 10: nn ≈ 8.0e7  (depleted as clouds form)
```
- ✅ **mean nn should DECREASE** over time
- ❌ If constant → re-initialization bug (should be fixed!)

#### **Cloud Formation (nc)**:
```
Expected pattern:
Step 1:  nc ≈ 10      (ncmin)
Step 2:  nc ≈ 250     (CCN activation starting)
Step 5:  nc ≈ 5.8e4   (clouds building)
Step 10: nc ≈ 2.0e5   (mature clouds)
```
- ✅ **mean nc should INCREASE** dramatically
- ✅ **max nc** can reach 1e5-1e6 #/kg in active clouds
- ❌ If starts at 5e7 → forced initialization bug (should be fixed!)

#### **Cloud Water (qc)**:
```
Expected pattern:
Step 1:  qc ≈ 0       (no clouds yet)
Step 2:  qc ≈ 0.01    (condensation starting)
Step 5:  qc ≈ 0.5     (clouds forming)
Step 10: qc ≈ 2.0     (mature clouds) g/kg
```
- ✅ **mean qc should INCREASE** as clouds form
- ✅ **max qc** typically 0.1-5 g/kg in active clouds
- ✅ Follows nc increase (more droplets → more water)

#### **Rain Formation (qr, nr)**:
```
Expected pattern (if conditions right):
Step 1-3:  qr ≈ 0, nr ≈ 0.01  (no rain yet)
Step 5:    qr ≈ 0.1, nr ≈ 100 (autoconversion starts)
Step 10:   qr ≈ 1.0, nr ≈ 1e4 (rain developing)
```
- ✅ Rain appears AFTER cloud water builds
- ✅ **qr and nr increase together**
- ✅ Autoconversion typically when qc > ~0.5 g/kg

#### **Ice/Snow/Graupel (qi, qs, qg)**:
```
Expected pattern (if cold enough):
- qi forms from deposition/freezing
- qs forms from ice aggregation
- qg forms from riming (ice + supercooled water)
```
- ✅ Ice species depend on temperature
- ✅ May stay zero if warm case
- ✅ Form gradually, slower than liquid

#### **Total Precipitation (Qp)**:
```
Qp = qr + qs + qg
```
- ✅ Should match sum of individual species
- ✅ Useful for tracking total precipitable water
- ✅ In ERF, this is computed automatically

---

## Comparison Checklist

### **Timestep 1 (Initialization)**:

| Variable | Check | Tolerance |
|----------|-------|-----------|
| qv | ERF ≈ WRF | < 0.1% |
| qc | ERF ≈ WRF | Exact if cold start |
| nn | ERF ≈ WRF | < 1% |
| nc | ERF ≈ WRF | Exact (~10) if cold start |
| nr | ERF ≈ WRF | Exact (~0.01) |

### **Timestep 5 (Evolution)**:

| Variable | Check | Tolerance |
|----------|-------|-----------|
| qv mean | ERF ≈ WRF | < 1% |
| qc mean | ERF ≈ WRF | < 5% |
| qc max | ERF ≈ WRF | < 10% |
| nn mean | Both decreasing | < 5% |
| nc mean | Both increasing | < 10% |
| nr mean | Both increasing | < 20% |

### **Timestep 10 (Mature State)**:

| Variable | Check | Tolerance |
|----------|-------|-----------|
| All q's | ERF ≈ WRF | < 5-10% |
| All n's | ERF ≈ WRF | < 10-20% |
| Trends | Same direction | Essential |

---

## Physical Sanity Checks

### **Typical Ranges** (for reference):

| Variable | Typical Range | Notes |
|----------|---------------|-------|
| **qv** | 1-15 g/kg | Depends on T, altitude |
| **qc** | 0-5 g/kg | Rarely > 3 g/kg |
| **qi** | 0-2 g/kg | Cold clouds only |
| **qr** | 0-10 g/kg | Heavy rain > 5 g/kg |
| **qs** | 0-5 g/kg | Snow/ice precip |
| **qg** | 0-5 g/kg | Hail/graupel |
| **nn** | 1e7-1e8 #/kg | Aerosol background |
| **nc** | 1e1-1e6 #/kg | 10-1M droplets/kg |
| **nr** | 1e-2-1e4 #/kg | Sparse to heavy rain |

### **Conservation Checks**:

1. **Total water conserved**:
   ```
   qv + qc + qi + qr + qs + qg = constant (minus precip loss)
   ```

2. **Number concentration consistency**:
   ```
   nc should decrease when: qc → qr (autoconversion)
   nc should increase when: nn → nc (CCN activation)
   ```

3. **Physical limits**:
   ```
   All q's ≥ 0 (no negative moisture!)
   All n's ≥ their minimum values (ncmin, nrmin, etc.)
   ```

---

## Troubleshooting

### **Problem: ERF and WRF qv differ significantly**

**Symptom**:
```
ERF: qv mean = 5.0 g/kg
WRF: qv mean = 8.0 g/kg
```

**Diagnosis**: Different initial conditions
- Check input files
- Verify same case is running
- Check if one includes moisture, other doesn't

---

### **Problem: ERF nc jumps to 5e7 immediately**

**Symptom**:
```
ERF Step 1: nc = 5.0e7 #/kg
WRF Step 1: nc = 1.0e1 #/kg
```

**Diagnosis**: Forced initialization bug (should be fixed!)
- Check that Fix #2 was applied correctly
- Verify ERF_InitWDM6.cpp has the simplified initialization

---

### **Problem: ERF nn stays constant, WRF nn decreases**

**Symptom**:
```
ERF: nn = 8.3e7 (all timesteps)
WRF: nn = 8.3e7 → 8.2e7 → 8.0e7
```

**Diagnosis**: Re-initialization bug (should be fixed!)
- Check that m_nn_initialized flag is working
- Verify it's only set once

---

### **Problem: Large differences in ice species**

**Symptom**:
```
ERF: qi = 0.5 g/kg
WRF: qi = 1.5 g/kg
```

**Diagnosis**: May be expected if using C++ GPU kernels
- C++ version doesn't have full ice physics yet
- Use Fortran bridge (`-DERF_ENABLE_WDM6_FORT=ON`) for exact match
- If already using Fortran, check temperature/initialization

---

### **Problem: Rain appears at different times**

**Symptom**:
```
ERF: qr > 0 at step 3
WRF: qr > 0 at step 5
```

**Diagnosis**: Could be numerical differences
- Check timestep sizes match
- Check autoconversion threshold
- Differences < 2 steps are acceptable
- Differences > 5 steps indicate problem

---

### **Problem: Statistics fluctuate wildly**

**Symptom**:
```
Step 5: nc = 1e5
Step 6: nc = 1e2
Step 7: nc = 1e6
```

**Diagnosis**: Possible instability
- Check timestep (may be too large)
- Check CFL condition
- Check for negative values being clipped
- May need smaller dt

---

## Example Comparison Script

Save as `compare_wdm6.py`:

```python
#!/usr/bin/env python3
import re
import sys

def parse_stats(filename):
    """Extract statistics from log file"""
    stats = {}
    with open(filename) as f:
        content = f.read()
        # Extract timestep blocks
        for match in re.finditer(r'Timestep\s+(\d+).*?========', content, re.DOTALL):
            step = int(match.group(1))
            block = match.group(0)
            stats[step] = {}
            # Extract each variable
            for var in ['qv', 'qc', 'qi', 'qr', 'qs', 'qg', 'Qp', 'nn', 'nc', 'nr']:
                if var in block:
                    # Extract mean value
                    m = re.search(rf'{var}.*?mean[=\s]+([\d.E+-]+)', block)
                    if m:
                        stats[step][var] = float(m.group(1))
    return stats

def compare(erf_file, wrf_file):
    """Compare ERF and WRF statistics"""
    erf = parse_stats(erf_file)
    wrf = parse_stats(wrf_file)
    
    print("Timestep | Variable | ERF Mean | WRF Mean | Diff % |")
    print("---------|----------|----------|----------|--------|")
    
    for step in sorted(erf.keys()):
        if step not in wrf:
            continue
        for var in ['qv', 'qc', 'nn', 'nc', 'nr', 'Qp']:
            if var in erf[step] and var in wrf[step]:
                erf_val = erf[step][var]
                wrf_val = wrf[step][var]
                if wrf_val != 0:
                    diff = 100 * abs(erf_val - wrf_val) / wrf_val
                    print(f"{step:8} | {var:8} | {erf_val:8.2e} | {wrf_val:8.2e} | {diff:6.2f} |")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: compare_wdm6.py <erf_log> <wrf_log>")
        sys.exit(1)
    compare(sys.argv[1], sys.argv[2])
```

**Usage**:
```bash
python3 compare_wdm6.py erf_wdm6_diag.log wrf_wdm6_diag.log
```

---

## Success Criteria

After comparing ERF and WRF outputs, you should see:

### **✅ Good Agreement**:
- qv matches within 1%
- qc, qr trends match (both increasing or both zero)
- nn decreases in both codes
- nc increases in both codes
- Absolute differences < 10% for active variables
- Same evolution patterns

### **✅ Expected Small Differences**:
- nc, nr can differ by ~10-20% (more variable)
- Ice species may differ if using C++ kernels
- Timing of events may differ by 1-2 timesteps
- Max values more variable than means

### **❌ Red Flags**:
- ERF nc starts at 5e7 (forced init bug)
- ERF nn constant (re-init bug)
- qv differs by > 5% (wrong initial conditions)
- Opposite trends (one increasing, one decreasing)
- Orders of magnitude difference

---

## Files Modified

### **ERF**:
- `Source/Microphysics/WDM6/ERF_UpdateWDM6.cpp` - Comprehensive statistics

### **WRF**:
- `/p/lustre1/wise14/wrf/compiling/WRF/WRF_v4.7.1_WDM6/phys/module_mp_wdm6.F` - Matching statistics

---

**Created**: 2026-07-31  
**By**: Claude Code (Sonnet 4)  
**Purpose**: Comprehensive hydrometeor comparison between ERF and WRF WDM6
