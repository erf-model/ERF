#!/usr/bin/env python3
"""
Phase 6 Two-Stream Radiation Time-Integration Consistency Check

Phase 6 audits the temporal semantics of radiation forcing in ERF's slow/fast
substep structure. This script validates that:

  1. Diagnostic CSV row count matches expected cadence (one per slow step)
  2. Heating rates remain finite and nonzero across all steps
  3. No unintended duplicate rows or missing rows due to timing errors
  4. Temporal consistency across a longer simulation (10 steps) vs Phase 5 (5 steps)

Unlike Phase 5 (which simply confirms wiring is active), Phase 6 focuses on
verifying that the temporal placement of radiation calls and source-term
injection produce expected diagnostic cadence without aliasing or repetition.
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

def check_timing_consistency():
    """Smoke-test Phase 6 temporal consistency of radiation forcing."""

    diag_file = "radiation_phase6_timing_diag.dat"
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

    # Test parameters (must match Phase6_TimeIntegration/inputs)
    S0 = 1361.0
    zenith_deg = 60.0
    cos_zenith = math.cos(math.radians(zenith_deg))
    expected_toa_flux = S0 * cos_zenith

    fixed_dt = 0.5
    stop_time = 5.0
    expected_min_steps = int(round(stop_time / fixed_dt))  # 10 steps

    print(f"\n{'='*70}")
    print("Phase 6 Two-Stream Radiation: Time-Integration Timing Check")
    print(f"{'='*70}")
    print(f"\nTest Parameters:")
    print(f"  Solar constant S0 = {S0:.2f} W/m^2")
    print(f"  Solar zenith angle = {zenith_deg:.1f}°")
    print(f"  cos(zenith) = {cos_zenith:.4f}")
    print(f"  Expected TOA flux = {expected_toa_flux:.4f} W/m^2")
    print(f"  Fixed dt = {fixed_dt} s, stop_time = {stop_time} s"
          f" (expect ~{expected_min_steps} diagnostic rows)")

    n_rows = len(data['step'])
    print(f"\nDiagnostic CSV rows found: {n_rows}")

    errors = []

    # ----
    # Check 1: Row count is as expected (not double, not missing)
    # ----
    # Phase 6 runs for 10 steps; expect ~10 rows if called once per step.
    # If 2x (20 rows) or 0.5x (5 rows), suspect temporal wiring error.
    max_expected = expected_min_steps + 2  # Allow small variation
    min_expected = expected_min_steps - 2
    if n_rows < min_expected or n_rows > max_expected:
        errors.append(
            f"Expected approximately {expected_min_steps} diagnostic rows "
            f"(one per slow step over {stop_time}s with dt={fixed_dt}s), "
            f"but found {n_rows}. This suggests advance_radiation() "
            f"may be called {n_rows / expected_min_steps:.1f}x per step. "
            f"Check for temporal wiring errors (Phase 6 D.12)."
        )
        print(f"  Row count check: {n_rows} rows (expected ~{expected_min_steps}) [FAIL]")
    else:
        print(f"  Row count check: {n_rows} rows ≈ {expected_min_steps} expected [PASS]")

    # ----
    # Check 2: All (step, time) pairs are unique (no duplicate rows)
    # ----
    step_time_pairs = list(zip(data['step'], data['time']))
    if len(step_time_pairs) != len(set(step_time_pairs)):
        errors.append(
            f"Duplicate (step, time) pairs detected in diagnostic CSV. "
            f"This may indicate advance_radiation() is called multiple "
            f"times per step, or diagnostics are being logged multiple times "
            f"per advance() call."
        )
        print(f"  Unique (step, time) pairs check [FAIL]")
    else:
        print(f"  Unique (step, time) pairs check [PASS]")

    # ----
    # Check 3: SW_TOA matches analytical value at every step
    # ----
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

    # ----
    # Check 4: Heating rate is finite and nonzero at every step
    # ----
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
            "nonzero SW+LW heating with sw_enabled=true, lw_enabled=true. "
            "This suggests heating rates are not being computed or "
            "injected into the source term."
        )
        print(f"  heating_rate_max nonzero check [FAIL]")
    else:
        print(f"  heating_rate_max nonzero check [PASS]")
        last_hr = data['heating_rate_max'][-1]
        mean_hr = sum(data['heating_rate_max']) / n_rows
        print(f"    (mean heating_rate_max = {mean_hr:.6e} K/s, "
              f"final = {last_hr:.6e} K/s)")

    # ----
    # Check 5: Heating rate is stable across the simulation
    # ----
    # Compute coefficient of variation (std dev / mean) to check for
    # unintended trends or instability
    if n_rows > 1 and nonzero_ok:
        import statistics
        mean_hr = statistics.mean(data['heating_rate_max'])
        stdev_hr = statistics.stdev(data['heating_rate_max'])
        cv_hr = stdev_hr / mean_hr
        print(f"  heating_rate_max stability (CV = {cv_hr:.4f})", end="")
        if cv_hr < 0.2:  # <20% variation
            print(" [PASS (stable)]")
        elif cv_hr < 0.5:  # <50% variation
            print(" [PASS (moderate variation)]")
        else:
            errors.append(
                f"heating_rate_max shows high variability (CV={cv_hr:.4f}). "
                f"Expected stable heating across all steps when using old-state "
                f"based forcing (Phase 6 R13). Check for transient startup "
                f"effects or unintended state-dependent radiation."
            )
            print(" [WARNING (high variation)]")

    # ----
    # Check 6: Other diagnostic columns remain finite
    # ----
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
        print("TEST PASSED - Phase 6 temporal consistency confirmed:")
        print("  - Diagnostic cadence matches expected (1 row per slow step)")
        print("  - Heating rates are finite, nonzero, and stable")
        print("  - No evidence of timing/wiring errors (Phase 6 D.12)")
        return True

if __name__ == "__main__":
    success = check_timing_consistency()
    sys.exit(0 if success else 1)
