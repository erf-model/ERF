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
    """Read the radiation diagnostic CSV and return a dict of column lists.

    The file is comma separated with a header line
    (step,time,call_site,SW_surface,SW_TOA,F_up_surface,F_down_toa,heating_rate_max,...),
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


def read_input_real(inputs_file, key, default):
    """Return the numeric value of `key` from an ERF inputs file, or `default`."""
    try:
        with open(inputs_file, 'r') as f:
            for line in f:
                line = line.split('#', 1)[0].strip()
                if not line or '=' not in line:
                    continue
                k, v = (s.strip() for s in line.split('=', 1))
                if k == key:
                    return float(v.strip('"'))
    except IOError:
        pass
    return default

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
    tau_per_layer = 0.003125  # Optical depth per layer [1/m]
    dz = 1024.0 / 64  # Layer thickness from domain setup: 64 layers over 10 km
    
    # Compute tau per unit height
    tau_per_m = tau_per_layer / dz  # Convert from per-layer to per-meter
    
    # Expected TOA flux (direct beam at top of atmosphere)
    expected_toa_flux = S0 * cos_zenith
    
    # Expected incident surface flux (after passing through all 64 layers)
    tau_total = tau_per_layer * 64
    expected_incident_flux = S0 * cos_zenith * math.exp(-tau_total / cos_zenith)

    # The SW_surface diagnostic is the flux absorbed by the surface, i.e. the
    # incident flux times (1 - albedo). The albedo comes from the inputs file
    # (erf.radiation.surface_albedo_sw, default 0.3 in RadChoice).
    surface_albedo = read_input_real("inputs", "erf.radiation.surface_albedo_sw", 0.3)
    expected_surface_flux = expected_incident_flux * (1.0 - surface_albedo)
    
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
    print(f"  Total optical depth (64 layers) = {tau_total:.4f}")
    print(f"\nExpected Fluxes:")
    print(f"  Expected TOA flux = {expected_toa_flux:.2f} W/m^2")
    print(f"  Surface albedo = {surface_albedo:.3f}")
    print(f"  Expected incident surface flux = {expected_incident_flux:.2f} W/m^2")
    print(f"  Expected absorbed surface flux = {expected_surface_flux:.2f} W/m^2")
    
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
