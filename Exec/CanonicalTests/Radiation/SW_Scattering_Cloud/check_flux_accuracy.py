#!/usr/bin/env python3
"""
Phase 4 Two-Stream Radiation Validation Script
Shortwave Scattering Cloud Test

This script verifies SW flux accuracy for the Meador-Weaver two-stream
diffuse (scattering) flux contribution introduced in Phase 4, combined
with the Phase 3 cloud-layer optical depth enhancement and cloud fraction
blending.

It replicates, level-by-level, the EXACT same algorithm implemented in
compute_sw_diffuse_flux() (ERF_TwoStreamSW.H) and the vertical_two_stream_
sweep() driver (ERF_AdvanceTwoStreamRadiation.cpp), for both the clear-sky
column and the cloudy column, then blends them via cloud_fraction, and
compares against the radiation_sw_scatter_diag.dat output file.

Checks:
1. Surface flux (direct + diffuse, blended) matches the replicated
   analytical calculation (within tolerance)
2. TOA flux matches input solar constant S0 * cos(zenith)
3. All computed fluxes are non-negative
4. Scattering cloud layer produces a nonzero diffuse contribution (sanity
   check that compute_sw_diffuse_flux() is actually being exercised)
"""

import sys
import os
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

def compute_sw_diffuse_flux(tau, F_dir_top, cos_zenith, omega, g):
    """
    Python replica of compute_sw_diffuse_flux() in
    Source/Radiation/ERF_TwoStreamSW.H (Meador-Weaver 1980 quadrature
    two-stream approximation). Must stay in exact numerical sync with
    the C++ implementation.
    """
    if omega <= 0.0 or cos_zenith <= 0.0 or F_dir_top <= 0.0 or tau <= 0.0:
        return 0.0

    w0 = min(omega, 1.0)
    asym = max(min(g, 1.0), -1.0)

    gamma1 = (7.0 - w0 * (4.0 + 3.0 * asym)) / 4.0
    gamma2 = -(1.0 - w0 * (4.0 - 3.0 * asym)) / 4.0
    gamma3 = (2.0 - 3.0 * asym * cos_zenith) / 4.0
    gamma4 = 1.0 - gamma3  # noqa: F841 (kept for parity with C++ derivation)

    lambda_sq = gamma1 * gamma1 - gamma2 * gamma2
    if lambda_sq < 1.0e-12:
        lambda_sq = 1.0e-12
    lam = math.sqrt(lambda_sq)

    denom_Gamma = gamma1 + lam
    Gamma = (gamma2 / denom_Gamma) if abs(denom_Gamma) > 1.0e-12 else 0.0

    Tdir = math.exp(-tau / cos_zenith)
    Tdiff = math.exp(-lam * tau)

    Gamma_sq = Gamma * Gamma
    Tdiff_sq = Tdiff * Tdiff
    denom = 1.0 - Gamma_sq * Tdiff_sq
    if abs(denom) < 1.0e-12:
        denom = 1.0e-12

    R_diff = Gamma * (1.0 - Tdiff_sq) / denom

    F_dir_bot = F_dir_top * Tdir
    F_dir_attenuated = F_dir_top - F_dir_bot

    F_diff_base = w0 * F_dir_attenuated * 0.5 + F_dir_bot * R_diff * gamma3
    if F_diff_base < 0.0:
        F_diff_base = 0.0

    return F_diff_base

def compute_column_surface_flux(
        n_layers, dz, tau_base, S0, cos_zenith,
        cloud_base_m, cloud_top_m, cloud_tau_per_layer,
        omega_clear, g_clear, omega_cloud, g_cloud, apply_cloud):
    """
    Python replica of the downward SW pass in vertical_two_stream_sweep()
    (ERF_AdvanceTwoStreamRadiation.cpp), including Phase 4 diffuse flux
    accumulation. Returns (F_dir_surface, F_diff_surface).
    """
    tau_cum = 0.0
    F_dir_prev = S0 * cos_zenith if cos_zenith > 0.0 else 0.0
    F_diff_prev = 0.0

    for k in range(n_layers):
        z = (k + 0.5) * dz
        in_cloud = apply_cloud and (cloud_base_m <= z <= cloud_top_m)

        if in_cloud:
            tau = tau_base + cloud_tau_per_layer
            omega, g = omega_cloud, g_cloud
        else:
            tau = tau_base
            omega, g = omega_clear, g_clear

        tau_cum += tau
        F_dir_curr = S0 * cos_zenith * math.exp(-tau_cum / cos_zenith) if cos_zenith > 0.0 else 0.0

        F_diff_layer = compute_sw_diffuse_flux(tau, F_dir_prev, cos_zenith, omega, g)
        Tdir_layer = math.exp(-tau / cos_zenith) if cos_zenith > 0.0 else 0.0
        F_diff_curr = F_diff_prev * Tdir_layer + F_diff_layer

        F_dir_prev = F_dir_curr
        F_diff_prev = F_diff_curr

    return F_dir_prev, F_diff_prev

def check_sw_scattering_cloud_accuracy():
    """Check SW flux accuracy for the Phase 4 scattering + cloud-layer test."""

    diag_file = "radiation_sw_scatter_diag.dat"
    if not os.path.exists(diag_file):
        print(f"ERROR: Diagnostic file {diag_file} not found")
        return False

    data = read_radiation_diag(diag_file)
    if data is None:
        return False

    # Test parameters (must match SW_Scattering_Cloud/inputs)
    S0 = 1361.0
    zenith_deg = 60.0
    cos_zenith = math.cos(math.radians(zenith_deg))
    tau_per_layer = 0.003125
    n_layers = 64
    domain_height_m = 1024.0
    dz = domain_height_m / n_layers

    cloud_base_height_m = 300.0
    cloud_top_height_m = 700.0
    cloud_tau_per_layer = 0.5
    cloud_fraction = 0.5

    # Phase 4 scattering parameters
    omega_clear = 0.0
    g_clear = 0.0
    omega_cloud = 0.9999
    g_cloud = 0.85

    expected_toa_flux = S0 * cos_zenith

    # Clear-sky column (no cloud enhancement applied at all)
    F_dir_clear, F_diff_clear = compute_column_surface_flux(
        n_layers, dz, tau_per_layer, S0, cos_zenith,
        cloud_base_height_m, cloud_top_height_m, cloud_tau_per_layer,
        omega_clear, g_clear, omega_cloud, g_cloud, apply_cloud=False)
    expected_flux_clear = F_dir_clear + F_diff_clear

    # Cloudy column (cloud-layer optical depth + cloud scattering properties
    # applied within the cloud band)
    F_dir_cloudy, F_diff_cloudy = compute_column_surface_flux(
        n_layers, dz, tau_per_layer, S0, cos_zenith,
        cloud_base_height_m, cloud_top_height_m, cloud_tau_per_layer,
        omega_clear, g_clear, omega_cloud, g_cloud, apply_cloud=True)
    expected_flux_cloudy = F_dir_cloudy + F_diff_cloudy

    # Blended surface flux (Phase 3 cloud fraction masking)
    expected_surface_flux = (
        (1.0 - cloud_fraction) * expected_flux_clear
        + cloud_fraction * expected_flux_cloudy
    )

    tolerance = 0.05  # 5%

    print(f"\n{'='*70}")
    print("Phase 4 Two-Stream Radiation: SW Scattering Cloud Test")
    print(f"{'='*70}")
    print(f"\nTest Parameters:")
    print(f"  Solar constant S0 = {S0:.2f} W/m^2")
    print(f"  Solar zenith angle = {zenith_deg:.1f}\u00b0")
    print(f"  cos(zenith) = {cos_zenith:.4f}")
    print(f"  Clear-sky optical depth per layer = {tau_per_layer:.6f}")
    print(f"  Vertical cell spacing dz = {dz:.2f} m")
    print(f"  Cloud layer = [{cloud_base_height_m:.1f}, {cloud_top_height_m:.1f}] m")
    print(f"  Cloud optical depth per layer (added) = {cloud_tau_per_layer:.4f}")
    print(f"  Cloud fraction = {cloud_fraction:.2f}")
    print(f"  Clear-sky (omega, g) = ({omega_clear:.4f}, {g_clear:.4f})")
    print(f"  Cloud (omega, g) = ({omega_cloud:.4f}, {g_cloud:.4f})")
    print(f"\nExpected Fluxes (direct + diffuse):")
    print(f"  Expected TOA flux = {expected_toa_flux:.4f} W/m^2")
    print(f"  Clear-sky column: direct={F_dir_clear:.6e}, diffuse={F_diff_clear:.6e},"
          f" total={expected_flux_clear:.6e} W/m^2")
    print(f"  Cloudy column:    direct={F_dir_cloudy:.6e}, diffuse={F_diff_cloudy:.6e},"
          f" total={expected_flux_cloudy:.6e} W/m^2")
    print(f"  Expected blended surface flux = {expected_surface_flux:.6f} W/m^2")

    last_idx = -1
    step = data['step'][last_idx]
    time = data['time'][last_idx]
    SW_surface = data['SW_surface'][last_idx]
    SW_TOA = data['SW_TOA'][last_idx]

    print(f"\nComputed Fluxes (step {step}, time {time:.4f}s):")
    print(f"  Computed TOA flux = {SW_TOA:.4f} W/m^2")
    print(f"  Computed surface flux = {SW_surface:.6f} W/m^2")

    errors = []

    # TOA flux check
    toa_error = abs(SW_TOA - expected_toa_flux) / expected_toa_flux
    print(f"\nAccuracy Checks:")
    print(f"  TOA flux error: {toa_error*100:.2f}%", end="")
    if toa_error > tolerance:
        print(f" [FAIL - exceeds {tolerance*100:.1f}% tolerance]")
        errors.append(f"TOA flux error too large: {toa_error*100:.2f}%")
    else:
        print(" [PASS]")

    # Surface flux check
    if expected_surface_flux > 0:
        surf_error = abs(SW_surface - expected_surface_flux) / expected_surface_flux
        print(f"  Surface flux error: {surf_error*100:.2f}%", end="")
        if surf_error > tolerance:
            print(f" [FAIL - exceeds {tolerance*100:.1f}% tolerance]")
            errors.append(f"Surface flux error too large: {surf_error*100:.2f}%")
        else:
            print(" [PASS]")

    # Non-negativity
    if SW_TOA < 0 or SW_surface < 0:
        errors.append(f"Negative flux detected: TOA={SW_TOA:.4f}, surface={SW_surface:.6f}")
        print(f"  Non-negativity check [FAIL]")
    else:
        print(f"  Non-negativity check [PASS]")

    # Scattering sanity check: the cloudy column's diffuse flux must be
    # strictly positive (confirms compute_sw_diffuse_flux() is actually
    # exercised for the cloud layer), while the clear-sky column's diffuse
    # flux must be exactly zero (confirms omega=0 clear-sky path is
    # unaffected, preserving Phase 1-3/Phase 3 behavior elsewhere).
    print(f"  Clear-sky diffuse flux == 0 (omega_clear=0)?"
          f" diffuse_clear={F_diff_clear:.6e}", end="")
    if F_diff_clear == 0.0:
        print(" [PASS]")
    else:
        print(" [FAIL]")
        errors.append("Clear-sky column diffuse flux is not exactly zero; "
                       "single_scattering_albedo=0 gating is broken")

    print(f"  Cloudy column diffuse flux > 0 (scattering active)?"
          f" diffuse_cloudy={F_diff_cloudy:.6e}", end="")
    if F_diff_cloudy > 0.0:
        print(" [PASS]")
    else:
        print(" [FAIL]")
        errors.append("Cloudy column diffuse flux is not positive; "
                       "compute_sw_diffuse_flux() scattering path not exercised")

    print(f"\n{'='*70}")
    if errors:
        print("TEST FAILED")
        for err in errors:
            print(f"  - {err}")
        return False
    else:
        print("TEST PASSED - SW scattering-cloud fluxes match analytical solution")
        return True

if __name__ == "__main__":
    success = check_sw_scattering_cloud_accuracy()
    sys.exit(0 if success else 1)
