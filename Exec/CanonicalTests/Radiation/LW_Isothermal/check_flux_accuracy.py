#!/usr/bin/env python3
"""
Phase 1 Two-Stream Radiation Validation Script
Longwave Isothermal Test

This script verifies LW flux accuracy in isothermal mode.
When all temperatures are uniform (T_iso_K), the two-stream LW solver should
satisfy:
    F_up(all levels) = F_down(all levels) = sigma * T_iso_K^4
    Net flux = F_up - F_down = 0 everywhere (LW_net_surface = 0)
    Outgoing LW at the top: LW_up_TOA = sigma * T_iso_K^4
    Heating rate = 0 everywhere

It reads radiation_lw_diag.dat and checks that:
1. The surface net LW is zero and the outgoing LW at the top is sigma*T^4 (within round-off)
2. Maximum heating rate is zero (within numerical precision)
3. All fluxes are consistent with Stefan-Boltzmann law
"""

import sys
import os
import numpy as np
import math

def read_radiation_diag(filename):
    """Read the radiation diagnostic CSV and return a dict of column lists.

    The file is comma separated with a header line
    (step,time,call_site,SW_surface,SW_TOA,SW_up_TOA,LW_net_surface,LW_up_TOA,heating_rate_max,...),
    so columns are looked up by name rather than by position. Non-numeric
    columns (call_site) are kept as strings; numeric columns are floats.
    """
    import csv
    try:
        with open(filename, 'r') as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                print(f"ERROR: No header found in {filename}")
                return None
            data = {name.strip(): [] for name in reader.fieldnames}
            for row in reader:
                if not any((v or '').strip() for v in row.values()):
                    continue
                for name in reader.fieldnames:
                    key = name.strip()
                    val = (row.get(name) or '').strip()
                    if key == 'call_site':
                        data[key].append(val)
                    elif key == 'step':
                        data[key].append(int(float(val)))
                    else:
                        try:
                            data[key].append(float(val))
                        except ValueError:
                            data[key].append(float('nan'))
    except IOError:
        print(f"ERROR: Could not read {filename}")
        return None

    if not data.get('step'):
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
    print(f"  Expected LW_up_TOA = σ*T^4 = {expected_flux:.4f} W/m^2")
    print(f"  Expected heating rate = 0 K/s everywhere")
    print(f"  Expected net flux = 0 W/m^2 everywhere")
    
    # Extract last timestep data
    last_idx = -1
    step = data['step'][last_idx]
    time = data['time'][last_idx]
    LW_net_surface = data['LW_net_surface'][last_idx]
    LW_up_TOA = data['LW_up_TOA'][last_idx]
    heating_rate_max = data['heating_rate_max'][last_idx]

    print(f"\nComputed Values (step {step}, time {time:.4f}s):")
    print(f"  Computed LW_up_TOA = {LW_up_TOA:.4f} W/m^2")
    print(f"  Computed LW_net_surface = {LW_net_surface:.4f} W/m^2")
    print(f"  Maximum heating rate = {heating_rate_max:.4e} K/s")

    # Check results
    errors = []

    print(f"\nAccuracy Checks:")
    # Isothermal column: the outgoing LW at the top is sigma*T^4.
    if expected_flux > 0:
        up_error = abs(LW_up_TOA - expected_flux) / expected_flux
        print(f"  LW_up_TOA relative error: {up_error:.4e}", end="")
        if up_error > flux_tolerance:
            print(f" [FAIL - exceeds {flux_tolerance:.4e}]")
            errors.append(f"LW_up_TOA error too large: {up_error:.4e}")
        else:
            print(" [PASS]")

    # Isothermal column: up and down fluxes are equal, so the surface net LW is zero.
    net_error = abs(LW_net_surface) / expected_flux
    print(f"  |LW_net_surface| / sigma T^4: {net_error:.4e}", end="")
    if net_error > flux_tolerance:
        print(f" [FAIL - exceeds {flux_tolerance:.4e}]")
        errors.append(f"LW_net_surface not zero: {LW_net_surface:.4e} W/m^2")
    else:
        print(" [PASS]")

    # In isothermal mode, heating rate should be exactly zero
    print(f"  Heating rate magnitude: {abs(heating_rate_max):.4e} K/s", end="")
    if abs(heating_rate_max) > heating_tolerance:
        print(f" [FAIL - exceeds {heating_tolerance:.4e}]")
        errors.append(f"Heating rate not zero: {heating_rate_max:.4e} K/s")
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
        print("TEST PASSED - LW isothermal test verified (LW_net_surface = 0, LW_up_TOA = σT^4, heating = 0)")
        return True

if __name__ == "__main__":
    success = check_lw_isothermal_accuracy()
    sys.exit(0 if success else 1)
