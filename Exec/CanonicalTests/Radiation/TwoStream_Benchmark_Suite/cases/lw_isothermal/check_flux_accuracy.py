#!/usr/bin/env python3
"""
Two-Stream Radiation Validation Script
Longwave Isothermal Column

The sounding prescribes theta(z) = T0 exp(g z / (c_p T0)) so that the
absolute temperature is T0 = 300 K at every level, and the surface is a
black body at the same T0. With a gray optical depth of 1 per layer over
64 layers the column is opaque, so the gray-gas solution gives

    LW_up_TOA      = sigma T0^4                (emission of the top layers)
    LW_net_surface = sigma T0^4 exp(-tau_col)  ~ 0  (F_down(0) -> sigma T0^4)

and every layer cools (cooling to space, strongest at the top), so
heating_rate_max is non-zero. This replaces the former isothermal_test
override, which forced these values instead of computing them.

It reads radiation_lw_diag.dat and checks that:
1. LW_up_TOA matches sigma*T0^4 within a small tolerance (the sounding is a
   piecewise-linear approximation of the exponential theta profile)
2. |LW_net_surface| is negligible compared to sigma*T0^4
3. heating_rate_max is finite and non-zero (the column radiates to space)
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
    
    # Test parameters (must match inputs file and sounding)
    T_iso_K = 300.0  # Isothermal temperature [K] (surface_temp_k and the sounding)
    sigma = 5.670374419e-8  # Stefan-Boltzmann constant [W/(m^2·K^4)]
    
    # Expected upwelling/downwelling flux (same in isothermal mode)
    expected_flux = compute_bb_radiation(T_iso_K, sigma)
    
    # Tolerance for numerical accuracy
    # In isothermal mode with exact arithmetic:
    #   - Upwelling and downwelling should be identical
    #   - Heating rates should be exactly zero
    # With finite precision, we allow ~1e-10 relative error
    # The 5-point sounding approximates the exponential theta profile to
    # ~0.05 K, i.e. ~0.1% in sigma T^4; allow 0.5%.
    flux_tolerance = 5e-3
    net_tolerance = 5e-3      # |LW_net_surface| / sigma T^4
    
    print(f"\n{'='*70}")
    print("Two-Stream Radiation: LW Isothermal Column Test")
    print(f"{'='*70}")
    print(f"\nTest Parameters:")
    print(f"  Isothermal temperature T_iso = {T_iso_K:.2f} K")
    print(f"  Stefan-Boltzmann constant σ = {sigma:.6e} W/(m^2·K^4)")
    print(f"\nAnalytical Solution:")
    print(f"  Expected LW_up_TOA = σ*T^4 = {expected_flux:.4f} W/m^2")
    print(f"  Expected LW_net_surface ≈ 0 W/m^2 (opaque column at the surface temperature)")
    print(f"  Expected heating_rate_max > 0 (cooling to space from the top layers)")
    
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
    if net_error > net_tolerance:
        print(f" [FAIL - exceeds {net_tolerance:.4e}]")
        errors.append(f"LW_net_surface not negligible: {LW_net_surface:.4e} W/m^2")
    else:
        print(" [PASS]")

    # The column radiates to space, so the heating (cooling) is finite and non-zero
    print(f"  Heating rate magnitude: {abs(heating_rate_max):.4e} K/s", end="")
    if not math.isfinite(heating_rate_max) or abs(heating_rate_max) <= 0.0:
        print(" [FAIL - expected a finite, non-zero cooling rate]")
        errors.append(f"heating_rate_max not finite/non-zero: {heating_rate_max}")
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
        print("TEST PASSED - isothermal column verified (LW_up_TOA = σT^4, LW_net_surface ≈ 0, cooling to space)")
        return True

if __name__ == "__main__":
    success = check_lw_isothermal_accuracy()
    sys.exit(0 if success else 1)
