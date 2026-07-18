#!/usr/bin/env python3
"""
Unit tests for the Albini (1983) firebrand spotting model.

Verifies lofting height formula, trajectory distance scaling, and
intensity threshold behaviour against Albini (1983) INT-309 reference values.

Reference values from:
  Albini, F.A. (1983). Potential spotting distance from wind-driven
    surface fires. USDA Forest Service Research Paper INT-309.
    Table 1: H_z values for I_B = 100, 500, 1000, 5000 kW/m.

Run: python3 test_albini_spotting.py
"""

import math
import sys


def albini_lofting_height(I_B_kW, I_B_min=100.0):
    """H_z = 12.2 * I_B^(1/3)  [Albini 1983 INT-309, eq. 1]"""
    if I_B_kW < I_B_min:
        return 0.0
    return 12.2 * (I_B_kW ** (1.0 / 3.0))


def albini_flight_time(H_z_m, v_t_ms=0.5):
    """t_f = H_z / v_t"""
    return H_z_m / v_t_ms if v_t_ms > 0 else 0.0


def albini_spot_distance_uniform_wind(I_B_kW, u_ms, v_t_ms=0.5):
    """
    Analytical spotting distance for uniform horizontal wind.
    d = u * t_f = u * H_z / v_t
    """
    H_z = albini_lofting_height(I_B_kW)
    t_f = albini_flight_time(H_z, v_t_ms)
    return u_ms * t_f


def test_lofting_height_threshold():
    """H_z = 0 below I_B_min threshold."""
    failures = []
    if albini_lofting_height(50.0, I_B_min=100.0) != 0.0:
        failures.append("I_B=50 kW/m should give H_z=0 (below threshold)")
    if albini_lofting_height(99.9, I_B_min=100.0) != 0.0:
        failures.append("I_B=99.9 kW/m should give H_z=0 (below threshold)")
    if albini_lofting_height(100.0, I_B_min=100.0) <= 0.0:
        failures.append("I_B=100 kW/m should give H_z > 0 (at threshold)")
    return failures


def test_lofting_height_formula():
    """
    Verify H_z = 12.2 * I_B^(1/3) against Albini (1983) INT-309 Table 1.
    Reference values (approximate):
      I_B=100 kW/m  -> H_z ~ 56 m
      I_B=500 kW/m  -> H_z ~ 96 m
      I_B=1000 kW/m -> H_z ~ 121 m
    """
    tol = 2.0  # 2 m tolerance on reference values
    cases = [(100.0, 56.0), (500.0, 96.0), (1000.0, 121.0)]
    failures = []
    for I_B, H_ref in cases:
        H_z = albini_lofting_height(I_B)
        if abs(H_z - H_ref) >= tol:
            failures.append(
                f"I_B={I_B} kW/m: H_z={H_z:.1f} m, expected ~{H_ref} m (tol={tol} m)"
            )
    return failures


def test_lofting_height_monotone():
    """H_z must increase monotonically with I_B."""
    I_values = [100.0, 200.0, 500.0, 1000.0, 5000.0]
    H_values = [albini_lofting_height(I) for I in I_values]
    failures = []
    for k in range(len(H_values) - 1):
        if H_values[k] >= H_values[k + 1]:
            failures.append(
                f"H_z not monotone: H({I_values[k]})={H_values[k]:.2f} >= "
                f"H({I_values[k+1]})={H_values[k+1]:.2f}"
            )
    return failures


def test_spot_distance_uniform_wind():
    """
    Spotting distance in uniform wind: d = u * H_z / v_t.
    For I_B=500 kW/m, u=5 m/s, v_t=0.5 m/s:
      H_z ~ 96 m, t_f = 96/0.5 = 192 s, d = 5 * 192 = 960 m
    """
    failures = []
    d = albini_spot_distance_uniform_wind(500.0, u_ms=5.0, v_t_ms=0.5)
    if not (800.0 < d < 1200.0):
        failures.append(f"Spot distance {d:.1f} m outside expected 800-1200 m range")
    return failures


def test_spot_distance_scales_with_wind():
    """Spotting distance scales linearly with wind speed."""
    failures = []
    d1 = albini_spot_distance_uniform_wind(500.0, u_ms=2.0)
    d2 = albini_spot_distance_uniform_wind(500.0, u_ms=4.0)
    ratio = d2 / d1
    if abs(ratio - 2.0) >= 0.01:
        failures.append(
            f"Distance ratio {ratio:.4f} should be 2.0 (linear in wind speed)"
        )
    return failures


def test_spot_distance_zero_below_threshold():
    """No spotting when I_B below threshold."""
    failures = []
    d = albini_spot_distance_uniform_wind(50.0, u_ms=10.0, v_t_ms=0.5)
    if d != 0.0:
        failures.append(f"Spot distance below threshold should be 0, got {d}")
    return failures


def test_terminal_velocity_effect():
    """Higher terminal velocity reduces flight time and spot distance."""
    failures = []
    d_slow = albini_spot_distance_uniform_wind(500.0, u_ms=5.0, v_t_ms=0.3)
    d_fast = albini_spot_distance_uniform_wind(500.0, u_ms=5.0, v_t_ms=1.0)
    if d_slow <= d_fast:
        failures.append(
            f"Slow fall ({d_slow:.1f} m) should give farther spotting than "
            f"fast fall ({d_fast:.1f} m)"
        )
    return failures


def main():
    print("=" * 70)
    print("Albini (1983) Spotting Model Unit Tests")
    print("Reference: Albini INT-309, USDA Forest Service, 1983")
    print("=" * 70)

    tests = [
        ("test_lofting_height_threshold",    test_lofting_height_threshold),
        ("test_lofting_height_formula",      test_lofting_height_formula),
        ("test_lofting_height_monotone",     test_lofting_height_monotone),
        ("test_spot_distance_uniform_wind",  test_spot_distance_uniform_wind),
        ("test_spot_distance_scales_with_wind", test_spot_distance_scales_with_wind),
        ("test_spot_distance_zero_below_threshold", test_spot_distance_zero_below_threshold),
        ("test_terminal_velocity_effect",    test_terminal_velocity_effect),
    ]

    failed = []
    for name, test_fn in tests:
        try:
            failures = test_fn()
            if failures:
                for msg in failures:
                    print(f"✗ {name} FAILED: {msg}")
                failed.append((name, failures[0]))
            else:
                print(f"✓ {name} PASSED")
        except Exception as e:
            print(f"✗ {name} ERROR: {e}")
            failed.append((name, str(e)))

    print("=" * 70)
    if not failed:
        print(f"All {len(tests)} tests PASSED")
        return 0
    else:
        print(f"{len(failed)} test(s) FAILED:")
        for name, error in failed:
            print(f"  - {name}: {error}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
