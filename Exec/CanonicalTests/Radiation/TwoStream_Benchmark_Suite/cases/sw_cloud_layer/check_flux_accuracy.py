#!/usr/bin/env python3
"""
Two-Stream Radiation Validation Script
Shortwave Cloud-Layer Test

Verifies the direct-beam SW surface flux for a column with a prescribed
cloud layer (tau_profile_type = cloud_layer) blended with a clear-sky
column through the cloud fraction:

    F_clear  = S0 cos(z) exp(-tau_clear / cos(z))
    F_cloudy = S0 cos(z) exp(-tau_cloudy / cos(z))
    F_sfc    = (1 - alpha) [ (1 - cf) F_clear + cf F_cloudy ]

where tau_cloudy adds cloud_tau_per_layer for every layer whose center lies
inside [cloud_base_height_m, cloud_top_height_m], alpha is the surface
albedo and cf the cloud fraction. No scattering is active in this test
(single_scattering_albedo = 0), so the diffuse contribution is exactly zero.

The parameters are read from the `inputs` file in the current directory so
the check stays consistent with the case configuration.
"""

import csv
import math
import os
import sys


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


def read_input_ints(inputs_file, key, default):
    """Return a list of ints for `key` (e.g. amr.n_cell), or `default`."""
    try:
        with open(inputs_file, 'r') as f:
            for line in f:
                line = line.split('#', 1)[0].strip()
                if not line or '=' not in line:
                    continue
                k, v = (s.strip() for s in line.split('=', 1))
                if k == key:
                    return [int(x) for x in v.split()]
    except IOError:
        pass
    return default


def read_input_reals(inputs_file, key, default):
    """Return a list of floats for `key` (e.g. geometry.prob_extent), or `default`."""
    try:
        with open(inputs_file, 'r') as f:
            for line in f:
                line = line.split('#', 1)[0].strip()
                if not line or '=' not in line:
                    continue
                k, v = (s.strip() for s in line.split('=', 1))
                if k == key:
                    return [float(x) for x in v.split()]
    except IOError:
        pass
    return default


def read_input_string(inputs_file, key, default):
    """Return the string value of `key` from an ERF inputs file, or `default`."""
    try:
        with open(inputs_file, 'r') as f:
            for line in f:
                line = line.split('#', 1)[0].strip()
                if not line or '=' not in line:
                    continue
                k, v = (s.strip() for s in line.split('=', 1))
                if k == key:
                    return v.strip('"')
    except IOError:
        pass
    return default


def read_radiation_diag(filename):
    """Read the radiation diagnostic CSV into a dict of column lists."""
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


def column_direct_flux(n_layers, dz, tau_base, S0, cos_zenith,
                       cloud_base_m, cloud_top_m, cloud_tau_per_layer, apply_cloud):
    """Direct-beam flux at the surface for a clear or cloudy column."""
    if cos_zenith <= 0.0:
        return 0.0
    tau_cum = 0.0
    for k in range(n_layers):
        z = (k + 0.5) * dz
        tau = tau_base
        if apply_cloud and cloud_base_m <= z <= cloud_top_m:
            tau += cloud_tau_per_layer
        tau_cum += tau
    return S0 * cos_zenith * math.exp(-tau_cum / cos_zenith)


def check_sw_cloud_layer_accuracy():
    inputs = "inputs"
    diag_file = read_input_string(inputs, "erf.radiation.diag_file", "radiation_sw_cloud_diag.dat")
    if not os.path.exists(diag_file):
        print(f"ERROR: Diagnostic file {diag_file} not found")
        return False

    data = read_radiation_diag(diag_file)
    if data is None:
        return False

    S0 = read_input_real(inputs, "erf.radiation.S0", 1361.0)
    zenith_deg = read_input_real(inputs, "erf.radiation.solar_zenith", 45.0)
    tau_per_layer = read_input_real(inputs, "erf.radiation.tau_per_layer", 0.05)
    cloud_base_m = read_input_real(inputs, "erf.radiation.cloud_base_height_m", 500.0)
    cloud_top_m = read_input_real(inputs, "erf.radiation.cloud_top_height_m", 1000.0)
    cloud_tau_per_layer = read_input_real(inputs, "erf.radiation.cloud_tau_per_layer", 0.5)
    cloud_fraction = read_input_real(inputs, "erf.radiation.cloud_fraction", 0.0)
    surface_albedo = read_input_real(inputs, "erf.radiation.surface_albedo_sw", 0.3)
    n_cell = read_input_ints(inputs, "amr.n_cell", [8, 8, 64])
    prob_extent = read_input_reals(inputs, "geometry.prob_extent", [3000.0, 3000.0, 1024.0])
    n_layers = n_cell[2]
    dz = prob_extent[2] / n_layers
    cos_zenith = math.cos(math.radians(zenith_deg))

    expected_toa_flux = S0 * cos_zenith
    F_clear = column_direct_flux(n_layers, dz, tau_per_layer, S0, cos_zenith,
                                 cloud_base_m, cloud_top_m, cloud_tau_per_layer, False)
    F_cloudy = column_direct_flux(n_layers, dz, tau_per_layer, S0, cos_zenith,
                                  cloud_base_m, cloud_top_m, cloud_tau_per_layer, True)
    expected_incident = (1.0 - cloud_fraction) * F_clear + cloud_fraction * F_cloudy
    expected_surface_flux = (1.0 - surface_albedo) * expected_incident

    tolerance = 0.05

    print(f"\n{'='*70}")
    print("Two-Stream Radiation: SW Cloud-Layer Test")
    print(f"{'='*70}")
    print("\nTest Parameters (from inputs):")
    print(f"  S0 = {S0:.2f} W/m^2, zenith = {zenith_deg:.1f} deg, cos(zenith) = {cos_zenith:.4f}")
    print(f"  Layers = {n_layers}, dz = {dz:.2f} m, tau_per_layer = {tau_per_layer:.6f}")
    print(f"  Cloud band = [{cloud_base_m:.1f}, {cloud_top_m:.1f}] m, "
          f"cloud_tau_per_layer = {cloud_tau_per_layer:.3f}, cloud_fraction = {cloud_fraction:.2f}")
    print(f"  Surface albedo = {surface_albedo:.3f}")
    print("\nExpected Fluxes:")
    print(f"  TOA flux = {expected_toa_flux:.2f} W/m^2")
    print(f"  Clear-column surface (incident) = {F_clear:.4f} W/m^2")
    print(f"  Cloudy-column surface (incident) = {F_cloudy:.4f} W/m^2")
    print(f"  Blended absorbed surface flux = {expected_surface_flux:.4f} W/m^2")

    SW_TOA = data['SW_TOA'][-1]
    SW_surface = data['SW_surface'][-1]
    print(f"\nComputed Fluxes (step {data['step'][-1]}, {data['call_site'][-1]}):")
    print(f"  TOA flux = {SW_TOA:.2f} W/m^2")
    print(f"  Surface flux = {SW_surface:.4f} W/m^2")

    errors = []
    toa_error = abs(SW_TOA - expected_toa_flux) / expected_toa_flux
    print("\nAccuracy Checks:")
    print(f"  TOA flux error: {toa_error*100:.2f}%", end="")
    if toa_error > tolerance:
        print(f" [FAIL - exceeds {tolerance*100:.1f}% tolerance]")
        errors.append(f"TOA flux error too large: {toa_error*100:.2f}%")
    else:
        print(" [PASS]")

    surf_error = abs(SW_surface - expected_surface_flux) / expected_surface_flux
    print(f"  Surface flux error: {surf_error*100:.2f}%", end="")
    if surf_error > tolerance:
        print(f" [FAIL - exceeds {tolerance*100:.1f}% tolerance]")
        errors.append(f"Surface flux error too large: {surf_error*100:.2f}%")
    else:
        print(" [PASS]")

    # The cloud must reduce the surface flux relative to the clear-sky value.
    clear_only = (1.0 - surface_albedo) * F_clear
    print(f"  Cloud attenuation (surface < clear-sky {clear_only:.4f})", end="")
    if cloud_fraction > 0.0 and cloud_tau_per_layer > 0.0 and not SW_surface < clear_only:
        print(" [FAIL]")
        errors.append("Surface flux not reduced by the cloud layer")
    else:
        print(" [PASS]")

    if any(v != v for v in data['SW_surface']) or any(v != v for v in data['heating_rate_max']):
        errors.append("NaN found in diagnostics")
        print("  Finite diagnostics check [FAIL]")
    else:
        print("  Finite diagnostics check [PASS]")

    print(f"\n{'='*70}")
    if errors:
        print("TEST FAILED")
        for err in errors:
            print(f"  - {err}")
        return False
    print("TEST PASSED - SW cloud-layer surface flux matches the analytical blend")
    return True


if __name__ == "__main__":
    sys.exit(0 if check_sw_cloud_layer_accuracy() else 1)
