# WDM6 Verification Tests

This directory contains test cases for verifying the WDM6 microphysics implementation.

## Test Cases

### 1. Bubble Test (Convective Cell)
**Purpose:** Basic functionality test - verify WDM6 runs and produces reasonable fields

**What to check:**
- WDM6 initializes without errors
- Number concentrations (nn, nc, nr) remain non-negative
- Cloud droplets form (nc > 0 where qc > 0)
- Rain drops form (nr > 0 where qr > 0)
- Precipitation accumulates at surface
- No NaN or Inf values

**Run:**
```bash
./erf3d.exe inputs_bubble_wdm6
```

### 2. CCN Sensitivity Test
**Purpose:** Verify CCN activation responds correctly to maritime vs continental conditions

**Test A: Maritime CCN (low concentration)**
```bash
./erf3d.exe inputs_bubble_wdm6_maritime
```
Expected: Lower nc (~100-300 cm^-3), larger droplets

**Test B: Continental CCN (high concentration)**  
```bash
./erf3d.exe inputs_bubble_wdm6_continental
```
Expected: Higher nc (~500-1000 cm^-3), smaller droplets

**Verification:**
- Maritime should have ~3-5x fewer droplets than continental
- Check cloud droplet effective radius: maritime > continental
- Total qc should be similar (same LWP)

### 3. WSM6 Comparison
**Purpose:** Compare WDM6 vs WSM6 behavior

**Run both:**
```bash
./erf3d.exe inputs_bubble_wsm6      # Single-moment reference
./erf3d.exe inputs_bubble_wdm6      # Double-moment
```

**Expected differences:**
- Precipitation timing: WDM6 may produce rain faster (autoconversion depends on nc)
- Cloud structure: WDM6 more realistic droplet spectrum
- Similar total precipitation over long integration
- Similar cloud ice fields (both single-moment for ice)

**Quantitative checks:**
- Total precip difference < 20%
- Cloud liquid water path (LWP) similar
- Peak updraft velocity similar

### 4. Number Concentration Bounds
**Purpose:** Verify physical bounds are maintained

**Check in plotfiles:**
- `nc >= ncmin = 1e1 m^-3` everywhere nc > 0
- `nr >= nrmin = 1e-2 m^-3` everywhere nr > 0  
- `nn >= 0` (can deplete to zero as CCN activate)
- `nc/qc` gives reasonable effective radius (5-50 μm for stratocumulus)
- `nr/qr` gives reasonable rain drop size

### 5. Mass Conservation
**Purpose:** Verify total water is conserved

**Check:**
```bash
# Compare initial and final total water
# qt = qv + qc + qi + qr + qs + qg (sum of all moisture)
# Should be conserved to roundoff
```

## Analysis Scripts

### plot_wdm6_fields.py
Plots WDM6-specific fields:
- Cloud droplet number concentration (nc)
- Rain drop number concentration (nr)
- Aerosol number concentration (nn)
- Cloud droplet effective radius
- Comparison of qc vs nc

### compare_wsm6_wdm6.py
Compares WSM6 and WDM6 output:
- Precipitation rates
- Cloud properties
- Thermodynamic profiles

### check_mass_conservation.py
Verifies total water conservation

## Expected Results

### Typical Values

**Stratocumulus cloud:**
- nc: 100-500 cm^-3 (maritime) to 300-1000 cm^-3 (continental)
- Cloud droplet radius: 8-15 μm
- LWP: 50-300 g/m^2

**Convective cloud:**
- nc: 50-300 cm^-3 at cloud base
- nr: 1-100 L^-1 in rain shaft
- Peak rain rate: 10-50 mm/hr

**Deep convection:**
- Strong updrafts → high supersaturation → high nc
- Graupel production similar to WSM6
- Faster autoconversion due to broader droplet spectrum

## Debugging

If WDM6 crashes or produces unphysical results:

1. **Check initialization:**
   - Is ccn0 set? (default 100e6 m^-3)
   - Are initial nc, nn, nr positive?

2. **Check for NaN:**
   ```bash
   # Enable FPE trapping in inputs file
   amrex.fpe_trap_invalid = 1
   ```

3. **Check microphysics debug output:**
   ```
   # In inputs file
   erf.microphysics_debug = 1
   ```

4. **Compare with WSM6:**
   - If WDM6 fails but WSM6 works, issue is in number concentration logic
   - Check CCN activation, droplet nucleation

5. **Check number/mass ratios:**
   - Unrealistic ratios indicate numerical issues
   - nc too high → timestep too small
   - nc/qc gives effective radius - should be 5-50 μm

## Performance Metrics

**Expected performance (relative to WSM6):**
- Fortran mode (ERF_USE_WDM6_FORT=ON): ~1.2-1.5x slower (CPU only)
- C++ GPU mode: ~1.1-1.3x slower (more variables to advect)
- Memory: +50% (3 additional 3D arrays)

## References

- WRF WDM6: `/p/lustre1/wise14/wrf/compiling/WRF/WRF_v4.7.1_clean/phys/module_mp_wdm6.F`
- Lim & Hong (2010): Development of an effective double-moment cloud microphysics scheme
- Morrison et al. (2009): For number concentration validation
