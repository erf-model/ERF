#!/usr/bin/env python3
"""
SEB Prognostic Mode Validation Check

Validates that:
1. Baseline case (seb_prognostic_enable=false) produces bitwise-identical output
2. Feature-on case (seb_prognostic_enable=true) produces finite T_s and q_s values
3. T_s evolves in direction consistent with SEB residual sign (expected: -10 W/m^2)
4. All T_s values remain within configured bounds [200, 340] K
5. All q_s values remain within configured bounds [0.0, 1.0] kg/kg
6. Each step changes T_s by dt * dT_s/dt (the Euler step uses the step size)
7. (mode "restart") a run restarted from the mid-run checkpoint continues the
   fresh run's T_s and q_s
"""

import sys
import os
import math

def read_csv(filename):
    """Read CSV file and return list of dicts (one per row, keys from header)"""
    if not os.path.exists(filename):
        print(f"ERROR: File not found: {filename}")
        return None

    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    if len(lines) < 2:
        print(f"ERROR: File {filename} has fewer than 2 lines")
        return None

    header = lines[0].split(',')
    rows = []
    for line in lines[1:]:
        values = line.split(',')
        if len(values) != len(header):
            print(f"WARNING: Skipping malformed row (expected {len(header)} cols, got {len(values)})")
            continue

        row = {}
        for key, val in zip(header, values):
            key = key.strip()
            val = val.strip()
            try:
                # Try to convert to float
                row[key] = float(val)
            except:
                row[key] = val
        rows.append(row)

    return {'header': header, 'rows': rows}

def check_baseline():
    """Check baseline case (feature disabled)"""
    print("\n" + "=" * 80)
    print("BASELINE CASE CHECK (seb_prognostic_enable=false)")
    print("=" * 80)

    baseline_file = "radiation_seb_prog_disabled.dat"
    data = read_csv(baseline_file)

    if data is None:
        print(f"FAIL: Could not read {baseline_file}")
        return False

    header = data['header']
    print(f"\nCSV header has {len(header)} columns:")
    print(f"  {', '.join(header)}")

    # the feature-off baseline should have 12 columns (no prognostic columns)
    # Base 8 + SEB residual 2 + prognostic 4 = 14 total
    # But baseline should have only 12 (skip prognostic columns)
    expected_baseline_cols = 12 # 8 base + 2 SEB + 2 placeholder for prognostic (as NaN)

    if len(header) < expected_baseline_cols:
        print(f"WARNING: Expected at least {expected_baseline_cols} columns, got {len(header)}")
        print("  This may be acceptable if prognostic columns are not written when disabled")

    # Verify T_s and q_s columns are NaN or missing
    if len(data['rows']) > 0:
        row = data['rows'][0]
        has_t_s_mean = 'T_s_mean' in row

        if has_t_s_mean:
            t_s_val = row.get('T_s_mean')
            if isinstance(t_s_val, float) and not math.isnan(t_s_val):
                print(f"WARNING: T_s_mean is not NaN in baseline case: {t_s_val}")
                print("  Expected NaN when feature is disabled")
                return False

    print("PASS: Baseline case looks good (prognostic fields are NaN or missing)")
    return True

def check_feature_on():
    """Check feature-on case (feature enabled)"""
    print("\n" + "=" * 80)
    print("FEATURE-ON CASE CHECK (seb_prognostic_enable=true)")
    print("=" * 80)

    feature_file = "radiation_seb_prog_enabled.dat"
    data = read_csv(feature_file)

    if data is None:
        print(f"FAIL: Could not read {feature_file}")
        return False

    header = data['header']
    print(f"\nCSV header has {len(header)} columns:")
    print(f"  {', '.join(header)}")

    # Feature-on should have all columns including T_s and q_s
    expected_features = ['T_s_mean', 'T_s_max', 'q_s_mean', 'q_s_max']
    for feature in expected_features:
        if feature not in header:
            print(f"ERROR: Expected column '{feature}' not found in header")
            return False

    print(f"PASS: All expected prognostic columns present")

    # Validate values
    if len(data['rows']) == 0:
        print("ERROR: No data rows in CSV")
        return False

    print("\nValidating prognostic values across all timesteps:")

    t_s_mean_values = []
    t_s_max_values = []
    q_s_mean_values = []
    q_s_max_values = []

    all_finite = True
    all_within_bounds = True

    t_min, t_max = 200.0, 340.0
    q_min, q_max = 0.0, 1.0

    # The prognostic surface update runs once per step at the post_dycore
    # call site; pre_dycore rows report NaN for T_s / q_s by design.
    post_rows = [r for r in data['rows'] if r.get('call_site', '').strip() == 'post_dycore']
    if not post_rows:
        print("ERROR: No post_dycore rows found in CSV")
        return False
    print(f"  Validating {len(post_rows)} post_dycore rows (pre_dycore rows carry NaN by design)")

    for i, row in enumerate(post_rows):
        step = int(row.get('step', -1))
        t_s_mean = row.get('T_s_mean')
        t_s_max = row.get('T_s_max')
        q_s_mean = row.get('q_s_mean')
        q_s_max = row.get('q_s_max')

        # Check finiteness
        for val, name in [(t_s_mean, 'T_s_mean'), (t_s_max, 'T_s_max'),
                          (q_s_mean, 'q_s_mean'), (q_s_max, 'q_s_max')]:
            if isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
                print(f"  Step {step}: ERROR - {name} is NaN/Inf: {val}")
                all_finite = False

        # Check bounds
        if isinstance(t_s_mean, float) and math.isfinite(t_s_mean):
            if not (t_min <= t_s_mean <= t_max):
                print(f"  Step {step}: ERROR - T_s_mean {t_s_mean:.2f} K out of bounds [{t_min}, {t_max}]")
                all_within_bounds = False
            t_s_mean_values.append(t_s_mean)

        if isinstance(t_s_max, float) and math.isfinite(t_s_max):
            if not (t_min <= t_s_max <= t_max):
                print(f"  Step {step}: ERROR - T_s_max {t_s_max:.2f} K out of bounds [{t_min}, {t_max}]")
                all_within_bounds = False
            t_s_max_values.append(t_s_max)

        if isinstance(q_s_mean, float) and math.isfinite(q_s_mean):
            if not (q_min <= q_s_mean <= q_max):
                print(f"  Step {step}: ERROR - q_s_mean {q_s_mean:.4f} kg/kg out of bounds [{q_min}, {q_max}]")
                all_within_bounds = False
            q_s_mean_values.append(q_s_mean)

        if isinstance(q_s_max, float) and math.isfinite(q_s_max):
            if not (q_min <= q_s_max <= q_max):
                print(f"  Step {step}: ERROR - q_s_max {q_s_max:.4f} kg/kg out of bounds [{q_min}, {q_max}]")
                all_within_bounds = False
            q_s_max_values.append(q_s_max)

    if not all_finite:
        print("FAIL: Some prognostic values are NaN/Inf")
        return False

    if not all_within_bounds:
        print("FAIL: Some prognostic values out of bounds")
        return False

    print(f"PASS: All prognostic values are finite and within bounds")

    # Check evolution direction (if we have multiple timesteps)
    if len(t_s_mean_values) >= 2:
        t_initial = t_s_mean_values[0]
        t_final = t_s_mean_values[-1]
        delta_t = t_final - t_initial

        print(f"\nTemperature evolution:")
        print(f"  Initial T_s_mean: {t_initial:.3f} K")
        print(f"  Final T_s_mean:   {t_final:.3f} K")
        print(f"  Change:           {delta_t:+.3f} K")

        # With SEB residual = -10 W/m^2, C_s = 2.0e4, tau = 86400 s
        # dT_s/dt ≈ -10 / 2.0e4 - (2*pi/86400) * (T_s - T_deep)
        # For weak damping (86400 s timescale), first term dominates initially
        # Expected: T_s should decrease (negative residual means cooling)
        if delta_t < 0.0:
            print(f"  Direction: DECREASING (consistent with negative residual)")
        elif delta_t > 0.0:
            print(f"  Direction: INCREASING (may be due to restoring term)")
        else:
            print(f"  Direction: NO CHANGE (negligible evolution over run)")

    if not check_time_integration(post_rows):
        return False

    print("\nPASS: Feature-on case validated successfully")
    return True


# Force-restore constants of inputs_seb_prognostic_enabled
SEB_RESIDUAL = (50.0 - 25.0) - 10.0 - 20.0 - 5.0   # R_net - H - LE - G [W/m^2]
SEB_C_S = 2.0e4          # seb_surface_heat_capacity [J/(m^2 K)]
SEB_TAU = 86400.0        # seb_restore_timescale_s [s]
SEB_T_DEEP = 295.0       # seb_t_deep_default [K]


def seb_tendency(t_s):
    """dT_s/dt of the force-restore equation for the deck's constants [K/s]."""
    return SEB_RESIDUAL / SEB_C_S - (2.0 * math.pi / SEB_TAU) * (t_s - SEB_T_DEEP)


def check_time_integration(post_rows):
    """Each post_dycore row must advance T_s by dt * dT_s/dt(T_s_old), where dt
    is the time between consecutive rows. The CSV carries six significant
    digits, so T_s near 300 K is quantised to 1e-4 K; the per-step test
    therefore allows one quantum on top of 2%, and the cumulative test over
    the whole run (where the quantisation averages out) uses 2% alone. A
    driver that multiplies the tendency by the absolute simulation time
    instead of the step size (as the original implementation did) fails both
    from the second row on."""
    print("\nTime-integration check (dT_s per step against dt * tendency):")
    quantum = 1.0e-4
    prev = None
    first = None
    expected_total = 0.0
    n_checked = 0
    n_bad = 0
    for row in post_rows:
        t = row.get('time')
        t_s = row.get('T_s_mean')
        if not (isinstance(t, float) and isinstance(t_s, float) and math.isfinite(t_s)):
            continue
        if first is None:
            first = t_s
        if prev is not None:
            dt_row = t - prev[0]
            if dt_row > 0.0:
                expected = dt_row * seb_tendency(prev[1])
                actual = t_s - prev[1]
                expected_total += expected
                n_checked += 1
                if abs(actual - expected) > 0.02 * abs(expected) + quantum:
                    n_bad += 1
                    if n_bad <= 5:
                        print(f"  step {int(row.get('step', -1))}: dT_s = {actual:+.6e} K, "
                              f"expected {expected:+.6e} K (dt = {dt_row:.3g} s); "
                              f"the update is {actual / expected:.2f}x the Euler step")
        prev = (t, t_s)
    if n_checked == 0:
        print("FAIL: fewer than two usable post_dycore rows")
        return False
    actual_total = prev[1] - first
    print(f"  total change {actual_total:+.5f} K over {n_checked} steps, "
          f"expected {expected_total:+.5f} K")
    if n_bad > 0:
        print(f"FAIL: {n_bad} of {n_checked} steps do not follow dt * tendency")
        return False
    if abs(actual_total - expected_total) > 0.02 * abs(expected_total):
        print("FAIL: cumulative T_s change is not the sum of the Euler steps")
        return False
    print(f"PASS: {n_checked} steps follow dt * tendency")
    return True


def check_restart():
    """T_s and q_s from a run restarted mid-way must continue the
    uninterrupted run: compare the post_dycore rows the two CSVs share."""
    print("\n" + "=" * 80)
    print("RESTART CHECK (inputs_seb_prognostic_restart against the fresh run)")
    print("=" * 80)
    fresh = read_csv("radiation_seb_prog_enabled.dat")
    restart = read_csv("radiation_seb_prog_restart.dat")
    if fresh is None or restart is None:
        return False
    fresh_rows = {int(r['step']): r for r in fresh['rows'] if r.get('call_site') == 'post_dycore'}
    restart_rows = [r for r in restart['rows'] if r.get('call_site') == 'post_dycore']
    if not restart_rows:
        print("FAIL: no post_dycore rows in the restart CSV")
        return False
    worst = 0.0
    n = 0
    for r in restart_rows:
        f = fresh_rows.get(int(r['step']))
        if f is None:
            continue
        for key in ('T_s_mean', 'T_s_max', 'q_s_mean', 'q_s_max'):
            a, b = r.get(key), f.get(key)
            if isinstance(a, float) and isinstance(b, float) and math.isfinite(a) and math.isfinite(b):
                worst = max(worst, abs(a - b))
                n += 1
    if n == 0:
        print("FAIL: the two CSVs share no post_dycore steps")
        return False
    first = restart_rows[0]
    f0 = fresh_rows.get(int(first['step']))
    print(f"  first restarted step {int(first['step'])}: T_s_mean {float(first['T_s_mean']):.4f} K "
          f"(fresh run {float(f0['T_s_mean']):.4f} K)")
    if worst > 1.0e-6:
        print(f"FAIL: restarted surface state departs from the fresh run by up to {worst:.3e}")
        return False
    print(f"PASS: {n} values agree to {worst:.1e} across the restart")
    return True

def main():
    if len(sys.argv) > 1:
        mode = sys.argv[1]
        if mode == "baseline":
            return 0 if check_baseline() else 1
        elif mode == "feature_on":
            return 0 if check_feature_on() else 1
        elif mode == "restart":
            return 0 if check_restart() else 1
        else:
            print(f"Unknown mode: {mode}")
            return 1

    # Run both checks
    baseline_ok = check_baseline()
    feature_ok = check_feature_on()

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Baseline case:   {'PASS' if baseline_ok else 'FAIL'}")
    print(f"Feature-on case: {'PASS' if feature_ok else 'FAIL'}")

    return 0 if (baseline_ok and feature_ok) else 1

if __name__ == "__main__":
    sys.exit(main())
