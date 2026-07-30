# Diagnosis: WDM6 Receiving No Moisture

## Problem Identified

Your diagnostic output shows:
```
WDM6::Advance() call #3600 (dt=1s)
  Number conc (#/kg): nn=[2009747687 to 2.020205033e+10], nc=[28627067.47 to 202020503.3], nr=[0 to 0]
  Max mixing ratios (g/kg): qc=0, qr=0, qi=0, qs=0, qg=0
```

**ALL mixing ratios are ZERO.** This is not a microphysics problem - WDM6 has no water to work with!

## Root Cause Options

Since CCN0 is correct (100e6 m^-3) and timestep is reasonable (1 second), the issue must be:

### Option 1: Initial conditions have no moisture
**Check:** Does your initial sounding or input file contain water vapor?

Typical supercell environment needs:
- Surface: qv ~ 14-18 g/kg (relative humidity ~ 70-80%)
- 500 mb: qv ~ 3-5 g/kg
- Without moisture, no clouds or precipitation can form!

**How to check:**
```bash
# Look at your initial conditions file or sounding
# Should have QVAPOR or mixing ratio specified
grep -i "qvapor\|mixing.*ratio\|moisture" your_input_file
```

### Option 2: Moisture exists but ERF isn't using it
**Check:** Is moisture actually being initialized in the ERF state?

The new diagnostics I added will print:
```
Copy_State_to_Micro diagnostic:
  Max qv in state = X g/kg
  Max qc in state = X g/kg
```

And in each WDM6 call:
```
  Max qv = X g/kg
```

**If these show 0.0 g/kg**, the problem is upstream of WDM6:
- Initial conditions not loaded properly
- Wrong variable names in input file
- Moisture model not enabled

**If these show > 0 g/kg**, there's a copy bug (unlikely based on the code I saw).

### Option 3: Saturation adjustment removing all moisture
**Check:** Is there a saturation adjustment or other physics removing moisture before WDM6 sees it?

Some models run saturation adjustment before microphysics, which could theoretically remove supersaturation. But it shouldn't remove ALL moisture.

## What to Do Next

### Step 1: Rebuild and rerun with new diagnostics

The changes I made will print:

1. **At first Copy_State_to_Micro call:**
   ```
   Copy_State_to_Micro diagnostic:
     Max qv in state = X g/kg
     Max qc in state = X g/kg
   ```

2. **Every WDM6::Advance() call now includes:**
   ```
     Max qv = X g/kg
   ```

Recompile and rerun. Look for these new lines in output.

### Step 2: Interpret the results

#### Case A: All qv values are 0
```
Copy_State_to_Micro diagnostic:
  Max qv in state = 0 g/kg
  Max qv = 0 g/kg
```

**Diagnosis:** Initial conditions have no moisture.

**Fix:** 
- Check your input file for moisture specification
- Verify QVAPOR is in initial condition file
- For idealized cases, ensure sounding has water vapor
- Typical supercell needs surface qv ~ 14-18 g/kg

#### Case B: State has moisture but WDM6 doesn't see it
```
Copy_State_to_Micro diagnostic:
  Max qv in state = 15.5 g/kg    <- moisture exists!
  Max qv = 0 g/kg                 <- but WDM6 doesn't see it!
```

**Diagnosis:** Bug in Copy_State_to_Micro (unlikely but possible).

**Fix:** 
- Check that moisture components (RhoQ1, RhoQ2, etc.) are correctly mapped
- Verify that `cons_in` is the right MultiFab

#### Case C: Moisture exists initially but disappears
```
Copy_State_to_Micro diagnostic:  (call #1)
  Max qv in state = 15.5 g/kg
  Max qv = 15.5 g/kg

WDM6::Advance() call #10:
  Max qv = 15.3 g/kg              <- moisture decreasing, OK
  
WDM6::Advance() call #100:
  Max qv = 0 g/kg                 <- moisture gone!
```

**Diagnosis:** Something is removing moisture between calls.

**Fix:**
- Check if other physics (radiation, PBL, etc.) are zeroing moisture
- Check for bugs in advection that might be causing numerical diffusion
- Verify boundary conditions aren't draining moisture

### Step 3: Check your input file

Please share or check:

1. **Initial condition type:**
   - Real data (WRF wrfinput file)?
   - Idealized (sounding-based)?
   - What format?

2. **Moisture variables in input:**
   ```bash
   # For WRF input files:
   ncdump -h wrfinput_d01 | grep -i "qvapor\|moisture"
   
   # For ERF input files:
   grep -i "moisture\|qvapor\|qv" your_input_file
   ```

3. **Physics options enabled:**
   ```bash
   grep -i "moisture\|microphysics\|mp_physics" your_input_file
   ```

4. **Is moisture model enabled in ERF?**
   - Check build configuration
   - Check input file for `erf.moisture_model = WDM6` or similar

## Quick Tests

### Test 1: Check if moisture is in PlotFiles
If you have earlier plotfiles (before t=3600s), check if they contain moisture:

```python
import yt
ds = yt.load("plt00001")  # early timestep
print(ds.field_list)      # should include ('boxlib', 'qv') or similar
print(ds.r['qv'].max())   # should be > 0
```

If early plotfiles show qv > 0 but later ones show qv = 0, moisture is being depleted.

### Test 2: Check domain moisture budget
```python
# For each plotfile, compute total domain moisture
import yt
for i in range(10):
    ds = yt.load(f"plt{i:05d}")
    qv_total = ds.r['qv'].sum()
    print(f"Step {i}: total qv = {qv_total}")
```

Should stay roughly constant (only decrease due to precipitation reaching surface).

### Test 3: Run simpler test case
Create minimal test with:
- Small domain (10x10x10 km)
- Uniform initial state with qv = 10 g/kg
- No winds initially
- Run for 10 timesteps

If WDM6 still sees qv = 0, the problem is in state initialization or copying.
If WDM6 sees qv > 0, the problem is specific to your supercell case setup.

## Most Likely Scenarios

Based on your output, I suspect one of these:

### Most likely: Initial conditions have no moisture
Many idealized test cases start dry and expect moisture to be added via surface fluxes or boundary conditions. For a supercell, you need:

1. **Moisture in initial sounding**
2. **Instability (CAPE > 1000 J/kg)**
3. **Shear (hodograph with veering)**

Check your initialization code/file for all three.

### Less likely: Wrong moisture index
If ERF is using a different component index for moisture, the copy would read the wrong array. But `RhoQ1_comp` is standard for qv.

### Unlikely: Saturation adjustment removed everything
This would require massive supersaturation removal, which doesn't make physical sense.

## What to Share

To help further, please provide:

1. **New diagnostic output** after rebuild:
   - The "Copy_State_to_Micro diagnostic" line
   - The "Max qv" line from WDM6::Advance
   
2. **How are you initializing the case?**
   - Real data input file?
   - Idealized sounding?
   - What's the source of initial conditions?

3. **Input file excerpt** showing:
   - Initial condition specification
   - Moisture/microphysics settings
   - Any relevant physics options

4. **First plotfile** (if available):
   - Can you check if qv > 0 at t=0?
   - Use yt or similar to examine

## Expected Outcome

After the diagnostic additions, you should see:

**Correct output:**
```
Copy_State_to_Micro diagnostic:
  Max qv in state = 15.5 g/kg
  Max qc in state = 0.001 g/kg

WDM6::Advance() call #1 (dt=1s)
  Max qv = 15.5 g/kg
  Number conc (#/kg): nn=[...], nc=[...], nr=[0 to 0]
  Max mixing ratios (g/kg): qc=0, qr=0, qi=0, qs=0, qg=0
  
WDM6::Advance() call #100 (dt=1s)
  Max qv = 15.2 g/kg          <- decreasing as clouds form
  Max mixing ratios (g/kg): qc=0.5, qr=0.01, ...
```

**If you still see qv=0**, we've confirmed it's an initialization problem, not a microphysics problem.
