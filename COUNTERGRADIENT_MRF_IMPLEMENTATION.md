# MRF PBL Scheme Countergradient Correction Implementation

## Overview

This implementation adds countergradient correction terms to the MRF (Medium-Range Forecast) PBL scheme in ERF. The countergradient correction is a key non-local turbulence transport term that captures updraft/downdraft correlations in unstable boundary layers.

## Technical Details

### 1. New EddyDiff Components

Added two new entries to the EddyDiff namespace enum in `Source/ERF_IndexDefines.H`:

- `CounterGradient_v`: Vertical countergradient correction term
- `CounterGradient_h`: Horizontal countergradient correction term (set to zero in current implementation)
- `NumDiffs`: Updated to account for new components (14 total)

### 2. New Configuration Parameter

Added `enable_mrf_countergradient` flag to `TurbChoice` struct in `Source/DataStructs/ERF_TurbStruct.H`:

- **Type**: `bool`
- **Default value**: `false` (maintains backward compatibility)
- **Usage**: Read from input via ParmParse with key `enable_mrf_countergradient`
- **Scope**: Per-level configurable (like other MRF parameters)

### 3. Physics Implementation

#### Formula
The countergradient correction term is computed as:

```
γ_c = b * (u_* * θ_*) / w_s
```

Where:
- `b` = pbl_mrf_const_b (default: 7.8)
- `u_*` = surface friction velocity (from surface layer model)
- `θ_*` = surface potential temperature scale (from surface layer model)
- `w_s` = u_* / φ_m (convective velocity scale)

#### Application
The countergradient term is only applied within the PBL (k < kpbl). In the free atmosphere and ghost cells, it is set to zero.

### 4. Modified Files

#### `Source/ERF_IndexDefines.H`
- Added `CounterGradient_v` and `CounterGradient_h` to EddyDiff enum
- Updated enum size from 11 to 14 components

#### `Source/DataStructs/ERF_TurbStruct.H`
- Added input parsing for `enable_mrf_countergradient` in MRF section
- Added member variable `bool enable_mrf_countergradient = false`
- Updated display() function to print new parameter

#### `Source/PBL/ERF_ComputeDiffusivityMRF.cpp`
- Added extraction of `enable_mrf_countergradient` flag before main computation loop
- Implemented countergradient calculation in PBL region
- Set countergradient terms to zero in free atmosphere
- Added countergradient terms to ghost cell initialization

## Usage

### Input File Example

To enable countergradient corrections with MRF PBL:

```
erf.pbl_type = "MRF"
erf.enable_mrf_countergradient = true
```

To disable (default):

```
erf.pbl_type = "MRF"
erf.enable_mrf_countergradient = false
```

### Output

When enabled, the countergradient correction values are stored in:
- `K_turb(i, j, k, EddyDiff::CounterGradient_v)` - vertical countergradient
- `K_turb(i, j, k, EddyDiff::CounterGradient_h)` - horizontal countergradient (currently zero)

These can be accessed in the diffusion routines and plotted for diagnostics.

## Backward Compatibility

The implementation maintains full backward compatibility:
- Default flag value is `false`
- When disabled, countergradient terms are set to zero
- Existing MRF tests and simulations are unaffected
- All other MRF functionality remains unchanged

## Future Enhancements

Potential improvements for future work:
1. Implement horizontal countergradient correction terms
2. Add entrainment parameterization at PBL top
3. Add cloud effects on stability
4. Implement entrainment flux terms following WRF
5. Add differential diffusion for moisture species

## References

- Hong and Pan (1996): "Nonlocal Boundary Layer Vertical Diffusion in a Medium-Range Forecast Model", MWR
- Hong et al. (2006): Further developments of MRF PBL scheme, MWR
- ERF Documentation: `/Docs/sphinx_doc/theory/PBLschemes.rst` - MRF PBL Model section
