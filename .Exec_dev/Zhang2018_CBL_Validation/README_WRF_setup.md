# WRF Setup for Zhang et al. (2018) SMS-3DTKE Validation

This directory contains a WRF namelist (`namelist.input.wrf`) configured to match the ERF SMS-3DTKE validation case.

## Overview

The case simulates a convective boundary layer (CBL) following Zhang et al. (2018) to validate the SMS-3DTKE turbulence closure scheme.

## Key Configuration

### Domain Setup
- **Horizontal**: 30 km × 30 km with 1 km resolution (32 × 32 grid points)
- **Vertical**: 2 km with ~80 m resolution (26 levels)
- **Boundary conditions**: Periodic in x and y, slip wall at top
- **Time step**: 5 seconds
- **Duration**: 4 hours (14400 seconds)

### Physics
- **Turbulence closure**: `km_opt = 5` (SMS-3DTKE - LES-type TKE scheme)
- **Subfilter-scale**: `sfs_opt = 2` (Subfilter TKE for SMS-3DTKE)
- **Diffusion**: `diff_opt = 2` (Full diffusion with 3D TKE)
- **PBL scheme**: `bl_pbl_physics = 0` (None - using LES closure instead)
- **Surface layer**: `sf_sfclay_physics = 1` (Revised MM5)
- **No radiation, microphysics, or cumulus schemes** (idealized case)

### SMS-3DTKE Parameters

The following parameters should match ERF configuration:

| Parameter | ERF Value | WRF Setting |
|-----------|-----------|-------------|
| C_k | 0.1 | `Ck` in SMS-3DTKE code |
| C_eps | 0.93 | `Ce` in SMS-3DTKE code |
| C_s | 0.18 | Smagorinsky constant |
| alpha_sgs | 0.25 | SGS fraction parameter |
| sigma_k | 0.5 | TKE Prandtl number |

**Note**: Some of these parameters may need to be hard-coded in WRF's SMS-3DTKE module or set via Registry modifications, as they are not all exposed as namelist parameters in standard WRF.

### Initial Condition

The initial sounding follows Zhang et al. (2018) Eq. (16):

```
z (m)      θ (K)         U (m/s)  V (m/s)
----------------------------------------------
0-925      300.0         10.0     0.0
925-1075   300 + 0.0536*(z-925)  10.0     0.0
>1075      308.05 + 0.003*(z-1075)  10.0   0.0
```

### Surface Forcing
- **Surface heat flux**: 0.24 K m/s (equivalent to ~300 W/m²)
- **Surface roughness**: z₀ = 0.16 m
- **No moisture flux**

### Large-Scale Forcing
- **Geostrophic wind**: U_g = 10 m/s, V_g = 0 m/s
- **Coriolis parameter**: f = 5×10⁻⁵ s⁻¹ (corresponds to ~45°N latitude)

## Setting Up WRF

### 1. Create Initial Condition

You'll need to create a custom `wrfinput_d01` file that matches the Zhang et al. (2018) sounding. Options:

**Option A: Use ideal.exe with custom module**
- Modify `dyn_em/module_initialize_ideal.F` to implement the Zhang sounding
- Compile and run `ideal.exe`

**Option B: Use WPS + custom met_em files**
- Create synthetic met_em files with the desired sounding
- Run `real.exe`

**Option C: Use Python/NCL to directly create wrfinput_d01**
- Start from a template `wrfinput` file
- Modify `T`, `U`, `V`, `P`, `PH` fields to match the sounding
- Set surface fields (TSK, HFX, etc.)

### 2. Set Surface Properties in wrfinput_d01

Key variables to set:
```python
# Surface roughness (all grid points)
wrfinput['ZNT'][0,:,:] = 0.16  # m

# Surface skin temperature (for heat flux calculation)
wrfinput['TSK'][0,:,:] = 300.0 + (0.24 / (Cd * U10))  # Adjust to get HFX ≈ 300 W/m²

# Land use category (to control surface scheme behavior)
wrfinput['LU_INDEX'][0,:,:] = 7  # Grassland or suitable category

# Soil moisture (if using land surface model)
# Set to avoid dry soil limiting heat flux
```

### 3. Verify SMS-3DTKE Implementation

Check that your WRF build includes SMS-3DTKE (km_opt = 5):

```bash
# In WRF/dyn_em directory, verify SMS-3DTKE code in module_diffusion.F
grep -i "sms.*3dtke\|km_opt.*5" dyn_em/module_diffusion.F

# Check Registry for km_opt and sfs_opt
grep "km_opt" Registry/Registry.EM_COMMON
grep "sfs_opt" Registry/Registry.EM_COMMON
```

If SMS-3DTKE (km_opt = 5) is not available in your WRF version, you'll need WRF 4.4+ or:
1. Obtain the SMS-3DTKE implementation patches
2. Apply to `dyn_em/module_diffusion.F` and related files
3. Update Registry entries if needed
4. Recompile WRF

### 4. Adjust Parameters (if needed)

Some SMS-3DTKE parameters may be hard-coded in `module_diffusion.F`. To match ERF exactly:

```fortran
! In dyn_em/module_diffusion.F, look for SMS-3DTKE parameters:
REAL, PARAMETER :: c_k = 0.1       ! Should match erf.Ck
REAL, PARAMETER :: c_eps = 0.93    ! Should match erf.Ce
REAL, PARAMETER :: c_s = 0.18      ! Should match erf.sms3dtke.C_s
REAL, PARAMETER :: alpha_sgs = 0.25  ! Should match erf.sms3dtke.alpha_sgs
REAL, PARAMETER :: prandtl = 0.5   ! Should match erf.sigma_k (or 1/sigma_k)
```

Note: Parameter names may vary by WRF version. Check Zhang et al. (2018) paper for exact definitions.

### 5. Run WRF

```bash
# Link namelist
ln -sf namelist.input.wrf namelist.input

# Run WRF
mpirun -np 4 ./wrf.exe
```

### 6. Monitor Output

Key diagnostics to compare with ERF:

- **Boundary layer height** (PBLH in WRF output)
- **Surface heat flux** (HFX)
- **Vertical profiles of**:
  - Potential temperature (T + perturbation)
  - Horizontal velocity (U, V)
  - TKE
  - Eddy diffusivity (K_m, K_h)
  - Momentum flux (u'w', v'w')
  - Heat flux (w'θ')

## Expected Results

Based on Zhang et al. (2018) and the ERF simulations:

- **PBL growth**: ~1000 m after 4 hours
- **Mixed layer**: Well-mixed θ profile up to PBL top
- **Entrainment zone**: Sharp gradient at PBL top
- **Surface flux**: Maintains ~0.24 K m/s throughout simulation
- **Wind profile**: Develops typical PBL structure with mixing

## Comparing ERF and WRF

Use the diagnostic output files from both models:

### ERF Output
- `mean_profiles.dat` - horizontally-averaged vertical profiles
- `flux_profiles.dat` - resolved and SGS flux profiles
- `sfs_profiles.dat` - SGS/subfilter-scale statistics
- `surf_hist.dat` - surface time series

### WRF Output
Extract from `wrfout*` files:
- Time-height cross sections of domain-averaged fields
- Surface time series from 2D fields
- TKE and mixing coefficients from PBL scheme output

### Key Metrics for Comparison

1. **PBL height evolution**: Should match within ~10%
2. **Temperature profile**: Mixed layer and entrainment zone structure
3. **TKE profiles**: Peak location and magnitude
4. **Surface fluxes**: Time series should be comparable
5. **Wind profiles**: Hodograph and jet structure

## Troubleshooting

### Heat flux not matching target
- Adjust TSK in wrfinput to control surface-air temperature difference
- Check ISFTCFLX and ISFFLX settings
- Verify surface layer scheme is active

### Simulation crashes
- Reduce time step if CFL violations occur
- Check for reasonable initial conditions (no NaNs, extreme gradients)
- Enable debug output in namelist

### Different PBL evolution
- Verify all SMS-3DTKE parameters match between codes
- Check vertical resolution matches (~80 m near surface)
- Confirm geostrophic forcing is applied correctly
- Verify Coriolis parameter (f = 5e-5 s⁻¹)

## References

Zhang, Y., W. I. Gustafson Jr., and J. D. Fast, 2018: Development and testing of a new scale-aware subgrid turbulent kinetic energy scheme. *Mon. Wea. Rev.*, **146**, 907–925, https://doi.org/10.1175/MWR-D-17-0220.1

## Contact

For questions about the ERF implementation, see the ERF repository.
For WRF-specific issues, consult WRF documentation or user forums.
