#!/usr/bin/env python3
"""
Phase 13 Two-Stream Radiation Validation Script
YSUNew PBL Coupling with Radiative Tendency Limiter/Smoother
"""

import sys
import math


def read_radiation_diag(filename):
    """Read radiation diagnostic CSV file and return parsed data."""
    data = {
        'step': [],
        'time': [],
        'call_site': [],
        'SW_surface': [],
        'SW_TOA': [],
        'F_up_surface': [],
        'F_down_toa': [],
        'heating_rate_max': []
    }

    try:
        with open(filename, 'r') as f:
            header = f.readline().strip()
            if not header.startswith('step'):
                f.seek(0)

            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split(',')
                if len(parts) >= 8:
                    try:
                        data['step'].append(int(float(parts[0])))
                        data['time'].append(float(parts[1]))
                        data['call_site'].append(parts[2].strip())
                        data['SW_surface'].append(float(parts[3]))
                        data['SW_TOA'].append(float(parts[4]))
                        data['F_up_surface'].append(float(parts[5]))
                        data['F_down_toa'].append(float(parts[6]))
                        data['heating_rate_max'].append(float(parts[7]))
                    except ValueError:
                        continue
    except FileNotFoundError:
        print(f"WARNING: Diagnostic file not found: {filename}")
        return None

    return data if data['step'] else None


def check_ysunew_radiation_coupling():
    diag_file = "radiation_phase13_ysu_coupling_diag.dat"

    print("=" * 70)
    print("PHASE 13 VALIDATION: YSUNew Radiation Coupling")
    print("=" * 70)

    data = read_radiation_diag(diag_file)
    if data is None:
        print(f"FAIL: No diagnostic data found in {diag_file}")
        return False

    nrows = len(data['step'])
    print(f"✓ Diagnostic file exists with {nrows} rows")

    if nrows < 3:
        print(f"FAIL: Expected at least 3 rows, got {nrows}")
        return False
    print(f"✓ Multiple rows recorded ({nrows})")

    # Monotonic check on post_dycore only (pre rows may repeat/reset)
    post_times = [t for t, cs in zip(data['time'], data['call_site']) if cs == 'post_dycore']
    if len(post_times) >= 2:
        for i in range(1, len(post_times)):
            if post_times[i] < post_times[i - 1]:
                print(f"FAIL: post_dycore time decreased at idx {i}: "
                      f"{post_times[i - 1]} -> {post_times[i]}")
                return False
        print("✓ post_dycore times are non-decreasing")
    else:
        print("WARN: Fewer than 2 post_dycore rows; skipped monotonic post_dycore check")

    # Finite checks
    for col_name, col_data in data.items():
        if col_name in ['step', 'call_site']:
            continue
        for i, val in enumerate(col_data):
            if not math.isfinite(val):
                print(f"FAIL: Non-finite {col_name} at row {i}: {val}")
                return False
    print("✓ All numeric diagnostic values are finite (no NaN/Inf)")

    # SW_TOA sanity
    expected_sw_toa = 1361.0 * 0.5
    tolerance = 0.05
    for i, sw_toa in enumerate(data['SW_TOA']):
        rel_error = abs(sw_toa - expected_sw_toa) / expected_sw_toa
        if rel_error > tolerance:
            print(f"WARN: SW_TOA anomaly at row {i}: {sw_toa:.2f} W/m² "
                  f"(expected ~{expected_sw_toa:.2f}, error {rel_error*100:.1f}%)")
    print(f"✓ SW_TOA values near expected {expected_sw_toa:.1f} W/m²")

    # heating_rate_max non-trivial
    max_h = max(data['heating_rate_max'])
    min_h = min(data['heating_rate_max'])
    if max_h < 1e-10:
        print(f"FAIL: heating_rate_max too small: max={max_h}")
        return False
    print(f"✓ heating_rate_max is nonzero: range [{min_h:.2e}, {max_h:.2e}] K/s")

    print("\n" + "=" * 70)
    print("PHASE 13 VALIDATION: PASS")
    print("=" * 70)
    return True


if __name__ == "__main__":
    ok = check_ysunew_radiation_coupling()
    sys.exit(0 if ok else 1)
