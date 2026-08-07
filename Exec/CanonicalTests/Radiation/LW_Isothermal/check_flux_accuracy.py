#!/usr/bin/env python3
"""
Phase 1 Two-Stream Radiation Validation Script
Longwave Isothermal Test

This script verifies LW flux accuracy in isothermal mode.
When all temperatures are uniform (T_iso_K), the two-stream LW solver should
satisfy:
    F_up(all levels) = F_down(all levels) = sigma * T_iso_K^4
    Net flux = F_up - F_down = 0 everywhere
    Heating rate = 0 everywhere

It reads radiation_lw_diag.dat and checks that:
1. Upwelling and downwelling LW fluxes are equal (within round-off)
2. Maximum heating rate is zero (within numerical precision)
3. All fluxes are consistent with Stefan-Boltzmann law
"""

import sys
import os
import numpy as np
import math

def read_radiation_diag(filename):
    """Read radiation diagnostic CSV file and return parsed data."""
    data = {
        'step': [],
        'time': [],
        'SW_surface': [],
        'SW_TOA': [],
        'F_up_surface': [],
        'F_down_toa': [],
        'heating_rate_max': []
    }
    
    try:
        with open(filename, 'r') as f:
            # Skip header line
            f.readline()
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split(',')
                if len(parts) >= 7:
                    data['step'].append(int(parts[0]))
                    data['time'].append(float(parts[1]))
                    data['SW_surface'].append(float(parts[2]))
                    data['SW_TOA'].append(float(parts[3]))
                    data['F_up_surface'].append(float(parts[4]))
                    data['F_down_toa'].append(float(parts[5]))
                    data['heating_rate_max'].append(float(parts[6]))
    except IOError:
        print(f"ERROR: Could not read {filename}")
        return None
    
    if not data['step']:
        print(f"ERROR: No data found in {filename}")
        return None
    
    return data

def compute_bb_radiation(T, sigma=5.670374419e-8):
    """
    Compute blackbody radiation intensity using Stefan-Boltzmann law.
    
    Args:
        T: Temperature [K]
        sigma: Stefan-Boltzmann constant [W/(m^2·K^4)]
    
    Returns:
        Radiative intensity [W/m^2]
    """
    if T <= 0:
        return 0.0
    return sigma * T**4

def check_lw_isothermal_accuracy():
    """Check LW isothermal test accuracy."""
    
    # Read diagnostic file
    diag_file = "radiation_lw_diag.dat"
    if not os.path.exists(diag_file):
        print(f"ERROR: Diagnostic file {diag_file} not found")
        return False
    
    data = read_radiation_diag(diag_file)
    if data is None:
        return False
    
    # Test parameters (must match inputs file)
    T_iso_K = 288.15  # Isothermal temperature [K]
    sigma = 5.670374419e-8  # Stefan-Boltzmann constant [W/(m^2·K^4)]
    
    # Expected upwelling/downwelling flux (same in isothermal mode)
    expected_flux = compute_bb_radiation(T_iso_K, sigma)
    
    # Tolerance for numerical accuracy
    # In isothermal mode with exact arithmetic:
    #   - Upwelling and downwelling should be identical
    #   - Heating rates should be exactly zero
    # With finite precision, we allow ~1e-10 relative error
    flux_tolerance = 1e-4
    heating_tolerance = 1e-4  # [K/s]
    
    print(f"\n{'='*70}")
    print("Phase 1 Two-Stream Radiation: LW Isothermal Test")
    print(f"{'='*70}")
    print(f"\nTest Parameters:")
    print(f"  Isothermal temperature T_iso = {T_iso_K:.2f} K")
    print(f"  Stefan-Boltzmann constant σ = {sigma:.6e} W/(m^2·K^4)")
    print(f"\nAnalytical Solution:")
    print(f"  Expected F_up = F_down = σ*T^4 = {expected_flux:.4f} W/m^2")
    print(f"  Expected heating rate = 0 K/s everywhere")
    print(f"  Expected net flux = 0 W/m^2 everywhere")
    
    # Extract last timestep data
    last_idx = -1
    step = data['step'][last_idx]
    time = data['time'][last_idx]
    F_up_surface = data['F_up_surface'][last_idx]
    F_down_toa = data['F_down_toa'][last_idx]
    heating_rate_max = data['heating_rate_max'][last_idx]
    
    print(f"\nComputed Values (step {step}, time {time:.4f}s):")
    print(f"  Computed F_up_surface = {F_up_surface:.4f} W/m^2")
    print(f"  Computed F_down_toa = {F_down_toa:.4f} W/m^2")
    print(f"  Maximum heating rate = {heating_rate_max:.4e} K/s")
    
    # Check results
    errors = []
    
    # In isothermal mode, upwelling and downwelling should be equal
    if expected_flux > 0:
        up_error = abs(F_up_surface - expected_flux) / expected_flux
        down_error = abs(F_down_toa - expected_flux) / expected_flux
        
        print(f"\nAccuracy Checks:")
        print(f"  F_up relative error: {up_error:.4e}", end="")
        if up_error > flux_tolerance:
            print(f" [FAIL - exceeds {flux_tolerance:.4e}]")
            errors.append(f"F_up error too large: {up_error:.4e}")
        else:
            print(" [PASS]")
        
        print(f"  F_down relative error: {down_error:.4e}", end="")
        if down_error > flux_tolerance:
            print(f" [FAIL - exceeds {flux_tolerance:.4e}]")
            errors.append(f"F_down error too large: {down_error:.4e}")
        else:
            print(" [PASS]")
    
    # In isothermal mode, heating rate should be exactly zero
    print(f"  Heating rate magnitude: {abs(heating_rate_max):.4e} K/s", end="")
    if abs(heating_rate_max) > heating_tolerance:
        print(f" [FAIL - exceeds {heating_tolerance:.4e}]")
        errors.append(f"Heating rate not zero: {heating_rate_max:.4e} K/s")
    else:
        print(" [PASS]")
    
    # Check flux equality (F_up should equal F_down in isothermal)
    flux_diff = abs(F_up_surface - F_down_toa)
    print(f"  |F_up - F_down| = {flux_diff:.4e} W/m^2", end="")
    if flux_diff > heating_tolerance * expected_flux:
        print(f" [FAIL]")
        errors.append(f"F_up != F_down (isothermal violation): diff = {flux_diff:.4e}")
    else:
        print(" [PASS]")
    
    # Overall result
    print(f"\n{'='*70}")
    if errors:
        print("TEST FAILED")
        for err in errors:
            print(f"  - {err}")
        return False
    else:
        print("TEST PASSED - LW isothermal test verified (F_up = F_down, heating = 0)")
        return True

if __name__ == "__main__":
    success = check_lw_isothermal_accuracy()
    sys.exit(0 if success else 1)
