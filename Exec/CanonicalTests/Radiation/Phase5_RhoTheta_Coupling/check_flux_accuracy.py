#!/usr/bin/env python3
"""
Phase 5 Two-Stream Radiation Validation Script
RhoTheta Coupling Smoke Test

Unlike the Phase 1-4 RegTest checkers (which independently re-derive the
expected flux/heating values via a Python replica of the C++ algorithm),
this script validates the Phase 5 WIRING itself -- i.e. that:

  1. compute_twostream_radiation_diagnostics() is now actually being
     called every step from advance_radiation() (Phase 5 Step 3). Before
     Phase 5, this function was never called anywhere in the codebase, so
     the diagnostic CSV would not exist / would not accumulate multiple
     rows under a real simulation run.
  2. qheating_rates gets populated with finite, physically reasonable
     per-level heating rates every step (Phase 5 Step 1), which is a
     necessary (though not sufficient, since this script only reads the
     domain-averaged/max diagnostics, not the full 3D MultiFab) condition
     for the RhoTheta source-term injection (Phase 5 Step 4) to have a
     real effect.
  3. No NaN/Inf appears in any diagnostic column across multiple steps,
     which would indicate a numerical breakdown introduced by the new
     per-level heating-rate computation (e.g. the LW net-flux-divergence
     calculation added in Phase 5 Step 1, which is new code exercised here
     for the first time with lw_enabled=true and isothermal_test=false
     together).

This is intentionally a "smoke test" rather than a flux-accuracy check:
Phase 5 does not change the underlying flux formulas validated by
SW_ClearSky_Analytical / LW_Isothermal / SW_Cloud_Layer /
SW_Scattering_Cloud; it only wires the existing (already-validated)
per-level heating-rate calculation into the simulation's source terms.
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

def check_rhotheta_coupling_smoke_test():
    """Smoke-test the Phase 5 RhoTheta coupling wiring."""

    diag_file = "radiation_phase5_coupling_diag.dat"
    if not os.path.exists(diag_file):
        print(f"ERROR: Diagnostic file {diag_file} not found")
        print("This likely means compute_twostream_radiation_diagnostics() "
              "was never called -- check that advance_radiation() wiring "
              "(Phase 5 Step 3) is present and erf.radiation_type is set "
              "correctly in the inputs file.")
        return False

    data = read_radiation_diag(diag_file)
    if data is None:
        return False

    # Test parameters (must match Phase5_RhoTheta_Coupling/inputs)
    S0 = 1361.0
    zenith_deg = 60.0
    cos_zenith = math.cos(math.radians(zenith_deg))
    expected_toa_flux = S0 * cos_zenith

    fixed_dt = 0.5
    stop_time = 2.5
    expected_min_steps = int(round(stop_time / fixed_dt))  # 5 steps

    print(f"\n{'='*70}")
    print("Phase 5 Two-Stream Radiation: RhoTheta Coupling Smoke Test")
    print(f"{'='*70}")
    print(f"\nTest Parameters:")
    print(f"  Solar constant S0 = {S0:.2f} W/m^2")
    print(f"  Solar zenith angle = {zenith_deg:.1f}\u00b0")
    print(f"  cos(zenith) = {cos_zenith:.4f}")
    print(f"  Expected TOA flux = {expected_toa_flux:.4f} W/m^2")
    print(f"  Fixed dt = {fixed_dt} s, stop_time = {stop_time} s"
          f" (expect >= {expected_min_steps} diagnostic rows)")

    n_rows = len(data['step'])
    print(f"\nDiagnostic CSV rows found: {n_rows}")

    errors = []

    # --------------------------------------------------------------------
    # Check 1: Multiple rows present -- confirms the driver is being
    # called repeatedly across the simulation (Phase 5 Step 3 wiring),
    # not just once (or never, prior to Phase 5).
    # --------------------------------------------------------------------
    if n_rows < expected_min_steps:
        errors.append(
            f"Expected at least {expected_min_steps} diagnostic rows "
            f"(one per timestep), found only {n_rows}. This suggests "
            f"compute_twostream_radiation_diagnostics() is not being "
            f"called every step -- check the advance_radiation() wiring."
        )
        print(f"  Row count check: {n_rows} < {expected_min_steps} [FAIL]")
    else:
        print(f"  Row count check: {n_rows} >= {expected_min_steps} [PASS]")

    # --------------------------------------------------------------------
    # Check 2: SW_TOA matches the analytical value at every step (sanity
    # that the TOA diagnostic, which does not depend on the new per-level
    # heating-rate machinery, remains correct across repeated calls).
    # --------------------------------------------------------------------
    toa_ok = True
    tolerance = 0.01  # 1%
    for i in range(n_rows):
        toa_error = abs(data['SW_TOA'][i] - expected_toa_flux) / expected_toa_flux
        if toa_error > tolerance:
            toa_ok = False
            errors.append(
                f"Step {data['step'][i]}: SW_TOA error {toa_error*100:.2f}% "
                f"exceeds {tolerance*100:.1f}% tolerance "
                f"(got {data['SW_TOA'][i]:.4f}, expected {expected_toa_flux:.4f})"
            )
    print(f"  SW_TOA accuracy check (all {n_rows} steps)"
          f" {'[PASS]' if toa_ok else '[FAIL]'}")

    # --------------------------------------------------------------------
    # Check 3: heating_rate_max is finite and nonzero at every step --
    # confirms qheating_rates is genuinely being populated with real
    # (non-garbage, non-zero) values by the Phase 5 per-level heating-rate
    # computation, both for SW (existing since Phase 1-4, now written
    # per-level instead of reduced-only) and LW (newly computed in Phase 5
    # Step 1 via compute_lw_heating_rate(), previously dead code).
    # --------------------------------------------------------------------
    finite_ok = True
    nonzero_ok = True
    for i in range(n_rows):
        hr = data['heating_rate_max'][i]
        if math.isnan(hr) or math.isinf(hr):
            finite_ok = False
            errors.append(f"Step {data['step'][i]}: heating_rate_max is "
                           f"NaN/Inf ({hr})")
        if hr == 0.0:
            nonzero_ok = False

    print(f"  heating_rate_max finite (no NaN/Inf) check"
          f" {'[PASS]' if finite_ok else '[FAIL]'}")
    if not nonzero_ok:
        errors.append(
            "heating_rate_max is exactly zero at every step; expected "
            "nonzero SW+LW heating with sw_enabled=true, lw_enabled=true, "
            "isothermal_test=false. This suggests qheating_rates is not "
            "being populated (Phase 5 Step 1 regression) or the RhoTheta "
            "coupling wiring (Phase 5 Steps 3-4) is not active."
        )
        print(f"  heating_rate_max nonzero (at least one step) check [FAIL]")
    else:
        print(f"  heating_rate_max nonzero (at least one step) check [PASS]")
        last_hr = data['heating_rate_max'][-1]
        print(f"    (heating_rate_max at final step = {last_hr:.6e} K/s)")

    # --------------------------------------------------------------------
    # Check 4: no other diagnostic column has NaN/Inf (SW_surface,
    # F_up_surface, F_down_toa) -- broader numerical sanity check.
    # --------------------------------------------------------------------
    other_cols_ok = True
    for i in range(n_rows):
        for col in ('SW_surface', 'F_up_surface', 'F_down_toa'):
            v = data[col][i]
            if math.isnan(v) or math.isinf(v):
                other_cols_ok = False
                errors.append(f"Step {data['step'][i]}: {col} is NaN/Inf ({v})")
    print(f"  Other diagnostic columns finite check"
          f" {'[PASS]' if other_cols_ok else '[FAIL]'}")

    print(f"\n{'='*70}")
    if errors:
        print("TEST FAILED")
        for err in errors:
            print(f"  - {err}")
        return False
    else:
        print("TEST PASSED - Phase 5 RhoTheta coupling wiring confirmed "
              "active and numerically stable across multiple timesteps")
        return True

if __name__ == "__main__":
    success = check_rhotheta_coupling_smoke_test()
    sys.exit(0 if success else 1)
