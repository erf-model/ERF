#!/usr/bin/env python3
"""
Phase 1 Two-Stream Radiation Validation Script
Shortwave Clear-Sky Analytical Test

This script verifies SW flux accuracy against the Beer-Lambert analytical solution:
    F_dir(z) = S0 * cos(zenith) * exp(-tau_cumulative(z) / cos(zenith))

It reads the radiation_sw_diag.dat output file and checks that:
1. Surface flux matches analytical prediction (within numerical precision)
2. TOA flux matches input solar constant S0 * cos(zenith)
3. All computed fluxes are non-negative
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

def compute_analytical_sw_flux(z, z_surface, z_toa, S0, cos_zenith, tau_per_layer):
    """
    Compute analytical SW direct-beam flux at height z using Beer-Lambert law.
    
    Args:
        z: Height above surface [m]
        z_surface: Surface height [m]
        z_toa: TOA height [m]
        S0: Solar constant at TOA [W/m^2]
        cos_zenith: Cosine of solar zenith angle
        tau_per_layer: Optical depth per unit height [1/m]
    
    Returns:
        Downwelling SW flux at height z [W/m^2]
    """
    if cos_zenith <= 0:
        return 0.0
    
    # Optical depth from TOA to this level
    tau_cumulative = tau_per_layer * (z_toa - z)
    return S0 * cos_zenith * math.exp(-tau_cumulative / cos_zenith)

def check_sw_flux_accuracy():
    """Check SW flux accuracy against analytical solution."""
    
    # Read diagnostic file
    diag_file = "radiation_sw_diag.dat"
    if not os.path.exists(diag_file):
        print(f"ERROR: Diagnostic file {diag_file} not found")
        return False
    
    data = read_radiation_diag(diag_file)
    if data is None:
        return False
    
    # Test parameters (must match inputs file)
    S0 = 1361.0  # Solar constant [W/m^2]
    zenith_deg = 60.0  # Solar zenith angle [degrees]
    cos_zenith = math.cos(math.radians(zenith_deg))
    tau_per_layer = 0.05  # Optical depth per layer [1/m]
    dz = 10000.0 / 20  # Layer thickness from domain setup: 20 layers over 10 km
    
    # Compute tau per unit height
    tau_per_m = tau_per_layer / dz  # Convert from per-layer to per-meter
    
    # Expected TOA flux (direct beam at top of atmosphere)
    expected_toa_flux = S0 * cos_zenith
    
    # Expected surface flux (after passing through all 20 layers)
    tau_total = tau_per_layer * 20
    expected_surface_flux = S0 * cos_zenith * math.exp(-tau_total / cos_zenith)
    
    # Tolerance for numerical accuracy (5%)
    tolerance = 0.05
    
    print(f"\n{'='*70}")
    print("Phase 1 Two-Stream Radiation: SW Clear-Sky Analytical Test")
    print(f"{'='*70}")
    print(f"\nTest Parameters:")
    print(f"  Solar constant S0 = {S0:.2f} W/m^2")
    print(f"  Solar zenith angle = {zenith_deg:.1f}°")
    print(f"  cos(zenith) = {cos_zenith:.4f}")
    print(f"  Optical depth per layer = {tau_per_layer:.4f}")
    print(f"  Total optical depth (20 layers) = {tau_total:.4f}")
    print(f"\nExpected Fluxes:")
    print(f"  Expected TOA flux = {expected_toa_flux:.2f} W/m^2")
    print(f"  Expected surface flux = {expected_surface_flux:.2f} W/m^2")
    
    # Extract last timestep data
    last_idx = -1
    step = data['step'][last_idx]
    time = data['time'][last_idx]
    SW_surface = data['SW_surface'][last_idx]
    SW_TOA = data['SW_TOA'][last_idx]
    
    print(f"\nComputed Fluxes (step {step}, time {time:.4f}s):")
    print(f"  Computed TOA flux = {SW_TOA:.2f} W/m^2")
    print(f"  Computed surface flux = {SW_surface:.2f} W/m^2")
    
    # Check results
    errors = []
    
    # Check TOA flux (should match S0*cos(zenith))
    toa_error = abs(SW_TOA - expected_toa_flux) / expected_toa_flux
    print(f"\nAccuracy Checks:")
    print(f"  TOA flux error: {toa_error*100:.2f}%", end="")
    if toa_error > tolerance:
        print(f" [FAIL - exceeds {tolerance*100:.1f}% tolerance]")
        errors.append(f"TOA flux error too large: {toa_error*100:.2f}%")
    else:
        print(" [PASS]")
    
    # Check surface flux (should match analytical value)
    if expected_surface_flux > 0:
        surf_error = abs(SW_surface - expected_surface_flux) / expected_surface_flux
        print(f"  Surface flux error: {surf_error*100:.2f}%", end="")
        if surf_error > tolerance:
            print(f" [FAIL - exceeds {tolerance*100:.1f}% tolerance]")
            errors.append(f"Surface flux error too large: {surf_error*100:.2f}%")
        else:
            print(" [PASS]")
    
    # Check non-negativity
    if SW_TOA < 0 or SW_surface < 0:
        errors.append(f"Negative flux detected: TOA={SW_TOA:.2f}, surface={SW_surface:.2f}")
        print(f"  Non-negativity check [FAIL]")
    else:
        print(f"  Non-negativity check [PASS]")
    
    # Overall result
    print(f"\n{'='*70}")
    if errors:
        print("TEST FAILED")
        for err in errors:
            print(f"  - {err}")
        return False
    else:
        print("TEST PASSED - SW clear-sky fluxes match analytical solution")
        return True

if __name__ == "__main__":
    success = check_sw_flux_accuracy()
    sys.exit(0 if success else 1)
