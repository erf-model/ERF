# Solution: nc is 680× Too High!

## The Problem

Your diagnostic output shows:
```
nc = 6.77e9 #/kg = 6770 cm^-3  (at ρ=1 kg/m³)
```

**This is 680× higher than WRF's minimum!**

WRF defines:
```fortran
real, parameter, private :: ncmin = 1.e1   ! 10 #/kg minimum
```

At ρ = 1 kg/m³, this is:
- ncmin = 10 #/kg = **0.01 cm^-3** 

Wait, that can't be right... let me check the units!

## Unit Analysis

Actually, I need to check WRF's units more carefully. In WRF:

```fortran
nn(i,k,j) = ccn0    ! Line ~223
```

Where ccn0 is passed in. Let me check what units ccn0 has in WRF...

Looking at the WRF Registry:
- `ccn0` is aerosol concentration
- In WRF namelist: `CCN_NOM = 100.E6` (100 million per m³)
- But in the code, it's used as `nn = ccn0` directly

**Key insight:** If `nn = ccn0` and ccn0 = 1e8 m^-3, then **nn is in units of #/m³, NOT #/kg!**

But in your output:
```
nn = 1.14e10 #/kg
```

This suggests ERF is converting: `nn = ccn0 / rho`

Let me check your initialization code...

## Your ERF Initialization (ERF_InitWDM6.cpp)

Line 91:
```cpp
nn(i,j,k) = ccn0_local / rho(i,j,k);  // Converting m^-3 to kg^-1
```

So ERF uses **#/kg** units, but WRF uses **#/m³** units!

## Converting Your nc to Physical Units

Your values:
```
nc = 6.77e9 #/kg
rho ~ 1.0 kg/m³ (at surface)

Physical nc = 6.77e9 #/kg × 1.0 kg/m³ = 6.77e9 #/m³
            = 6770 #/cm³
```

**This is still way too high!**

Typical values:
- Clean maritime: 50-100 cm^-3 = 5e7 - 1e8 m^-3
- Continental: 200-500 cm^-3 = 2e8 - 5e8 m^-3
- Polluted: 500-1500 cm^-3 = 5e8 - 1.5e9 m^-3
- **Your value: 6770 cm^-3 = 6.77e9 m^-3** ← Extremely polluted!

## Root Cause in Your Code

In `ERF_InitWDM6.cpp` line 102:
```cpp
Real activation_fraction = Real(0.15);  // 15% activation
Real nc_init = nn(i,j,k) * activation_fraction;
```

With:
- nn = 1.14e10 #/kg
- activation_fraction = 0.15
- **nc = 1.7e9 #/kg** (from activation alone)

But your output shows nc = 6.77e9 #/kg, which is **4× higher** than even 15% activation!

**Something is adding to nc multiple times**, or the cap isn't working.

## WRF's Approach

WRF **does not initialize nc from nn** in the microphysics. Instead:

1. **Initial condition:** nc = 0 (no clouds)
2. **When clouds form:** nc increases through:
   - Nucleation processes (not in your simple init code)
   - Heterogeneous freezing that releases droplets
   - Or nc stays small and WRF uses diagnostic relationships

WRF's philosophy: **Let nc be small initially, let microphysics build it up naturally**

## The Fix

###  Option 1: Match WRF - Start with nc = 0 or very small

```cpp
if (nc(i,j,k) < Real(1.e1)) {
    if (qc(i,j,k) > Real(1.e-9)) {
        // Cloud exists - set small initial nc
        nc(i,j,k) = Real(1.e7) / rho(i,j,k);  // 10 cm^-3 = 1e7 m^-3
    } else {
        // No cloud - keep nc very small
        nc(i,j,k) = Real(1.e1);  // Matches WRF's ncmin
    }
}
```

### Option 2: Reduce activation fraction drastically

```cpp
Real activation_fraction = Real(0.001);  // 0.1% activation
Real nc_init = nn(i,j,k) * activation_fraction;
```

This would give:
- nc = 1.14e10 × 0.001 = 1.14e7 #/kg = 11 cm^-3 ← **Perfect!**

### Option 3: Set absolute cap

```cpp
Real nc_max = Real(5.0e8) / rho(i,j,k);  // 500 cm^-3 maximum
nc(i,j,k) = amrex::min(nc_init, nc_max);
```

## Recommended Solution

**Use Option 1 - Start very small like WRF:**

```cpp
// WDM6: Initialize cloud droplet number concentration
// Following WRF: start with very small nc, let microphysics build it up
if (nc(i,j,k) < Real(1.e1)) {
    if (qc(i,j,k) > Real(1.e-9)) {
        // Cloud already exists - set minimal realistic nc
        // 10 cm^-3 = 1e7 m^-3 / rho
        nc(i,j,k) = Real(1.e7) / rho(i,j,k);
    } else {
        // No cloud yet - use WRF's ncmin
        nc(i,j,k) = Real(1.e1);  // 10 #/kg = 0.01 cm^-3 at rho=1
    }
}
```

## Why This Will Help

With nc ~ 10-50 cm^-3 instead of 6770 cm^-3:

**Current (nc = 6770 cm^-3):**
- Mean diameter: D ~ 10-12 μm (tiny!)
- Autoconversion: Very slow (threshold not reached)
- Result: qr stays small, Qp ~ 1 g/kg

**After fix (nc = 10-50 cm^-3):**
- Mean diameter: D ~ 25-30 μm (normal!)
- Autoconversion: Active (larger drops collide efficiently)
- Result: qr should increase to 3-10 g/kg, Qp ~ 5-15 g/kg

## Mystery to Solve

Your nc = 6.77e9 #/kg is **4× higher** than 15% × nn.

Possible causes:
1. Initialization is being called multiple times
2. Something else is adding to nc
3. Background activation (line 115) is too high

Check if initialization is called more than once, or if there's other code modifying nc.

