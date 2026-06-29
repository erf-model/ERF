# MRF PBL Model Enhancements Testing

This directory contains test cases for the MRF (Medium Range Forecast) PBL model enhancements implemented in ERF. These enhancements align the MRF model more closely with WRF's advanced parameterizations while maintaining backward compatibility.

## Test Cases

### 1. Full_Enhancements
Tests the MRF model with all three enhancements enabled:
- Iterative thermal excess correction
- Cloud-aware moisture diffusivity (IMVDIF)
- Countergradient term bounding

This represents the "Full WRF-Style Configuration" for maximum realism in challenging conditions.

**Usage:**
```bash
./erf inputs_full_enhancements
```

### 2. Conservative_Enhancement
Tests the MRF model with only iterative thermal excess correction enabled. This conservative configuration improves strong convection without other changes.

**Usage:**
```bash
./erf inputs_conservative
```

### 3. CloudAware_Layers
Tests the MRF model with cloud-aware moisture diffusivity enabled. This configuration is suitable for stratocumulus-topped or fog/stratus simulations.

**Usage:**
```bash
./erf inputs_cloud_aware
```

### 4. Baseline_MRF
Tests the baseline MRF model without any enhancements. This serves as a reference for backward compatibility and neutral/stable boundary layer simulations.

**Usage:**
```bash
./erf inputs_baseline
```

## Enhancement Details

### Iterative Thermal Excess Correction
- **Parameter:** `pbl.enable_mrf_iterative_thermal_excess`
- **Sub-parameters:** `pbl.mrf_thermal_excess_iterations`
- **Physical basis:** Refines surface virtual potential temperature estimate through multiple iterations
- **Typical impact:** 1-5% additional computation

### Cloud-Aware Moisture Diffusivity (IMVDIF)
- **Parameter:** `pbl.enable_mrf_cloudy_layers`
- **Sub-parameters:** `pbl.mrf_cloud_diffusivity_factor` (default: 0.8)
- **Physical basis:** Reduces moisture diffusivity in stratocumulus layers
- **Typical impact:** <1% overhead

### Countergradient Term Bounding
- **Parameter:** `pbl.enable_mrf_countergradient_bounds`
- **Sub-parameters:** `pbl.mrf_countergradient_max_theta`, `pbl.mrf_countergradient_max_q`
- **Physical basis:** Prevents unrealistic amplification of countergradient fluxes
- **Typical impact:** <1% overhead

## Running Tests

All test cases use the same problem setup (uniform ABL simulation) with different PBL configurations. To run a test:

```bash
# Build ERF
cd Exec
make -j

# Run specific test
cd CanonicalTests/ABL/MRF_Enhancements
../../../erf inputs_full_enhancements
```

## Output and Validation

Each test produces plotfiles at intervals specified in the input files. The test cases are designed to run stably and demonstrate that all enhancements can be combined without causing numerical instability.

## References

- Hong et al. (1996): "Nonlocal Boundary Layer Vertical Diffusion in a Medium-Range Forecast Model"
- Hong et al. (2006): "A new vertical diffusion package with an explicit treatment of entrainment processes"
- WRF Model Documentation: PBL Schemes
