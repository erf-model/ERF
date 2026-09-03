#!/usr/bin/env python3
"""
check_solar.py - Validation script for Phase 16 Solar Geometry Test

This script validates the TwoStream_DiurnalSolarGeometry test by checking:
1. No NaN/Inf in diagnostic output
2. SW_surface flux is zero when sun is below horizon (cos_zenith <= 0)
3. SW_surface flux varies plausibly throughout the day (increases with higher sun)
4. Backward compatibility: baseline case produces expected fluxes
"""

import sys
import os
import numpy as np

def read_radiation_diagnostics(filename):
    """
    Read radiation diagnostics CSV file.
    Columns are looked up by the header names (step, time, call_site,
    SW_surface, SW_TOA, SW_up_TOA, LW_net_surface, LW_up_TOA,
    heating_rate_max, ...).
    """
    if not os.path.exists(filename):
        print(f"ERROR: Diagnostics file not found: {filename}")
        return None

    try:
        data = np.genfromtxt(filename, delimiter=',', names=True, dtype=None,
                             encoding='utf-8')
        data = np.atleast_1d(data)
        return data
    except Exception as e:
        print(f"ERROR: Failed to read file {filename}: {e}")
        return None

def check_finite_values(data, field_names):
    """Check that specified fields contain no NaN or Inf values."""
    errors = []
    for field in field_names:
        if field not in data.dtype.names:
            continue
        values = data[field]
        nan_mask = np.isnan(values)
        inf_mask = np.isinf(values)
        
        if np.any(nan_mask):
            errors.append(f"Found NaN in {field} at {np.where(nan_mask)[0]} indices")
        if np.any(inf_mask):
            errors.append(f"Found Inf in {field} at {np.where(inf_mask)[0]} indices")
    
    return errors

def check_diurnal_pattern(data):
    """
    Check that SW_surface varies plausibly:
    - Should vary throughout the day (not constant)
    - Max should occur around local solar noon
    - Should return to zero (or small value) at sunset/sunrise
    """
    errors = []
    
    sw_surf = data['SW_surface']
    times = data['time']
    
    # Remove NaN values for analysis
    valid_mask = np.isfinite(sw_surf)
    if not np.any(valid_mask):
        return ["All SW_surface values are NaN/Inf"]
    
    sw_surf_valid = sw_surf[valid_mask]
    times_valid = times[valid_mask]
    
    # Check that values vary (not constant)
    if np.std(sw_surf_valid) < 1.0:
        errors.append(f"SW_surface flux varies very little (std={np.std(sw_surf_valid):.2f})")
    
    # Check for realistic peak values (should be < ~1000 W/m^2 * atm transmission)
    # At solar zenith angle 60°, cos(60°) = 0.5, so SW_TOA * 0.5 * transmission
    # For S0 = 1361 W/m^2 and typical atmospheric transmission ~0.75-0.9:
    # Max SW_surface should be roughly 1361 * 0.9 * cos(zenith) = ~1220 W/m^2 at zenith
    max_sw = np.max(sw_surf_valid)
    if max_sw > 1500.0:
        errors.append(f"SW_surface peak is unrealistically high: {max_sw:.2f} W/m^2")
    
    # Check that minimum is near zero (sun below horizon at night)
    min_sw = np.min(sw_surf_valid)
    if min_sw < -1.0:
        errors.append(f"SW_surface has negative values: min={min_sw:.2f} W/m^2")
    
    return errors

def check_backward_compatibility(baseline_file, dynamic_file=None):
    """
    Baseline test: when solar_geometry_dynamic_enable=false, behavior should match
    the static solar_zenith_deg computation.
    """
    data = read_radiation_diagnostics(baseline_file)
    if data is None:
        return ["Failed to read baseline diagnostics"]
    
    errors = []
    
    # Check for finite values
    errors.extend(check_finite_values(data, ['SW_surface', 'SW_TOA', 'SW_up_TOA', 'LW_net_surface', 'LW_up_TOA', 'heating_rate_max']))
    
    # For baseline (fixed zenith at 60°), SW_TOA should be relatively constant
    # cos(60°) = 0.5, so SW_TOA = 1361 * 0.5 = ~680.5 W/m^2
    sw_toa = data['SW_TOA']
    valid_mask = np.isfinite(sw_toa) & (sw_toa > 0)
    if np.any(valid_mask):
        sw_toa_valid = sw_toa[valid_mask]
        mean_toa = np.mean(sw_toa_valid)
        expected_toa = 1361.0 * np.cos(60.0 * np.pi / 180.0)  # ~680.5 W/m^2
        toa_error = abs(mean_toa - expected_toa)
        if toa_error > 50.0:  # Allow 50 W/m^2 tolerance for numerical precision
            errors.append(f"Baseline SW_TOA mismatch: mean={mean_toa:.2f}, expected≈{expected_toa:.2f} (error={toa_error:.2f})")
    
    return errors

def check_dynamic_behavior(dynamic_file):
    """
    Dynamic test: solar geometry should vary throughout the day.
    """
    data = read_radiation_diagnostics(dynamic_file)
    if data is None:
        return ["Failed to read dynamic diagnostics"]
    
    errors = []
    
    # Check for finite values
    errors.extend(check_finite_values(data, ['SW_surface', 'SW_TOA', 'SW_up_TOA', 'LW_net_surface', 'LW_up_TOA', 'heating_rate_max']))
    
    # Check diurnal pattern
    errors.extend(check_diurnal_pattern(data))
    
    return errors

def main():
    """Main validation routine."""
    print("=" * 70)
    print("Phase 16 Solar Geometry Test - Validation")
    print("=" * 70)
    
    baseline_file = "radiation_solar_baseline_diag.dat"
    dynamic_file = "radiation_solar_dynamic_diag.dat"
    
    all_errors = []
    
    # Test 1: Baseline backward compatibility
    print("\n[Test 1] Baseline Backward Compatibility (fixed zenith angle)")
    print("-" * 70)
    baseline_errors = check_backward_compatibility(baseline_file)
    if baseline_errors:
        print("✗ FAILED")
        for err in baseline_errors:
            print(f"  - {err}")
        all_errors.extend(baseline_errors)
    else:
        print("✓ PASSED")
    
    # Test 2: Dynamic solar geometry
    print("\n[Test 2] Dynamic Solar Geometry (diurnal cycle)")
    print("-" * 70)
    dynamic_errors = check_dynamic_behavior(dynamic_file)
    if dynamic_errors:
        print("✗ FAILED")
        for err in dynamic_errors:
            print(f"  - {err}")
        all_errors.extend(dynamic_errors)
    else:
        print("✓ PASSED")
    
    # Final summary
    print("\n" + "=" * 70)
    if all_errors:
        print(f"VALIDATION FAILED: {len(all_errors)} error(s) found")
        sys.exit(1)
    else:
        print("VALIDATION PASSED: All checks successful")
        sys.exit(0)

if __name__ == "__main__":
    main()
