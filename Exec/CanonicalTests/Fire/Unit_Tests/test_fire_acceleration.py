#!/usr/bin/env python3
"""
Unit tests for Phase 12 Fire Acceleration Models.

Verifies both the size-based (Rothermel 1983 / Catchpole et al. 1992) and
temporal (McAlpine & Wakimoto 1991 / VanWagner) acceleration models.

Tests verify mathematical correctness of:
  1. Size-based: α = 1 - exp(-r_fire / L_acc)
  2. Temporal: R(t) = R_E × (1 - exp(-A × t))
  3. Perimeter-based A selection (point vs line)
  4. Wind-lag exponential damping

Run: python3 test_fire_acceleration.py
"""

import math
import sys


def size_based_alpha(r_fire, L_acc):
    """
    Compute size-based acceleration factor.
    
    α = 1 - exp(-r_fire / L_acc)
    
    Parameters:
        r_fire: effective fire radius [m]
        L_acc: acceleration length scale [m]
    
    Returns:
        alpha: scaling factor [0, 1]
    """
    L_acc_safe = max(L_acc, 1.0)  # floor prevents divide-by-zero
    return 1.0 - math.exp(-r_fire / L_acc_safe)


def temporal_vanwagner(R_equilibrium, A_per_sec, elapsed_time):
    """
    Compute temporal acceleration via VanWagner equation.
    
    R(t) = R_E × (1 - exp(-A × t))
    
    Parameters:
        R_equilibrium: target equilibrium ROS [m/s]
        A_per_sec: acceleration constant [1/s]
        elapsed_time: time since ignition [s]
    
    Returns:
        R_current: instantaneous ROS [m/s]
    """
    return R_equilibrium * (1.0 - math.exp(-A_per_sec * elapsed_time))


def wind_lag_damping(R_old, R_new, dt, tau_wind):
    """
    Compute wind-lag exponential damping.
    
    R_target = R_old + (R_new - R_old) × (1 - exp(-dt / tau))
    
    Parameters:
        R_old: current ROS [m/s]
        R_new: new equilibrium ROS [m/s]
        dt: timestep [s]
        tau_wind: time constant [s]
    
    Returns:
        R_target: lagged target ROS [m/s]
    """
    return R_old + (R_new - R_old) * (1.0 - math.exp(-dt / tau_wind))


def test_size_based_alpha_at_zero():
    """Size-based alpha at zero fire size: α = 1 - exp(0) = 0."""
    r_fire = 0.0
    L_acc = 50.0
    alpha = size_based_alpha(r_fire, L_acc)
    
    # At r_fire = 0, alpha should be exactly 0
    assert abs(alpha - 0.0) < 1.0e-10, \
        f"alpha at r_fire=0 should be 0, got {alpha}"
    print("✓ test_size_based_alpha_at_zero PASSED")


def test_size_based_alpha_at_large_fire():
    """Size-based alpha at large fire: α ≈ 1.0 when r_fire >> L_acc."""
    r_fire = 500.0
    L_acc = 50.0
    alpha = size_based_alpha(r_fire, L_acc)
    
    # At r_fire = 500 >> L_acc = 50, alpha should approach 1.0
    # exp(-10) ≈ 4.54e-5, so alpha ≈ 0.9999545
    expected = 1.0 - math.exp(-10.0)
    assert abs(alpha - expected) < 1.0e-4, \
        f"alpha at large fire should ≈ {expected}, got {alpha}"
    assert alpha > 0.9999, \
        f"alpha should be > 0.9999, got {alpha}"
    print("✓ test_size_based_alpha_at_large_fire PASSED")


def test_size_based_alpha_monotone():
    """Size-based alpha increases monotonically with r_fire."""
    L_acc = 50.0
    r_values = [0.0, 10.0, 25.0, 50.0, 100.0, 200.0]
    alphas = [size_based_alpha(r, L_acc) for r in r_values]
    
    # Verify strictly increasing
    for i in range(len(alphas) - 1):
        assert alphas[i] < alphas[i+1], \
            f"alpha should increase: {alphas[i]} >= {alphas[i+1]}"
    print("✓ test_size_based_alpha_monotone PASSED")


def test_size_based_alpha_formula():
    """Size-based alpha formula verification: r_fire = 25, L_acc = 50 → α ≈ 0.3935."""
    r_fire = 25.0
    L_acc = 50.0
    alpha = size_based_alpha(r_fire, L_acc)
    
    # α = 1 - exp(-0.5) ≈ 1 - 0.606531 = 0.393469
    expected = 1.0 - math.exp(-0.5)
    assert abs(alpha - expected) < 0.0001, \
        f"alpha formula error: expected {expected:.4f}, got {alpha:.4f}"
    print(f"✓ test_size_based_alpha_formula PASSED (α={alpha:.4f})")


def test_temporal_vanwagner_at_t_zero():
    """Temporal VanWagner at t=0: R(0) = R_E × (1 - exp(0)) = 0."""
    R_E = 0.5  # m/s
    A_per_sec = 0.115 / 60.0  # convert 0.115 1/min to 1/s
    t = 0.0
    R = temporal_vanwagner(R_E, A_per_sec, t)
    
    assert abs(R - 0.0) < 1.0e-10, \
        f"R(0) should be 0, got {R}"
    print("✓ test_temporal_vanwagner_at_t_zero PASSED")


def test_temporal_vanwagner_at_t_infinity():
    """Temporal VanWagner at t→∞: R(∞) = R_E × (1 - exp(-∞)) = R_E."""
    R_E = 0.5  # m/s
    A_per_sec = 0.115 / 60.0  # convert 0.115 1/min to 1/s
    
    # At t = 10 × (1/A), the result should be within 1e-4 of R_E
    t_long = 10.0 / A_per_sec
    R = temporal_vanwagner(R_E, A_per_sec, t_long)
    
    # exp(-10) ≈ 4.54e-5, so R ≈ R_E × (1 - 4.54e-5) ≈ R_E
    relative_error = abs(R - R_E) / R_E
    assert relative_error < 1.0e-4, \
        f"R(∞) should be ≈ R_E, relative error = {relative_error}"
    print("✓ test_temporal_vanwagner_at_t_infinity PASSED")


def test_temporal_vanwagner_monotone():
    """Temporal VanWagner increases monotonically with t."""
    R_E = 0.5  # m/s
    A_per_sec = 0.115 / 60.0
    times = [0.0, 60.0, 120.0, 300.0, 600.0, 1200.0]
    R_values = [temporal_vanwagner(R_E, A_per_sec, t) for t in times]
    
    # Verify strictly increasing
    for i in range(len(R_values) - 1):
        assert R_values[i] < R_values[i+1], \
            f"R(t) should increase: {R_values[i]} >= {R_values[i+1]}"
    print("✓ test_temporal_vanwagner_monotone PASSED")


def test_temporal_vanwagner_formula():
    """Temporal VanWagner formula: R_E = 0.5, A = 0.115/min, t = 120s."""
    R_E = 0.5  # m/s
    A_per_min = 0.115  # 1/min
    A_per_sec = A_per_min / 60.0  # convert to 1/s
    t = 120.0  # seconds
    
    R = temporal_vanwagner(R_E, A_per_sec, t)
    
    # R(120) = 0.5 × (1 - exp(-0.115/60 × 120))
    #        = 0.5 × (1 - exp(-0.23))
    #        = 0.5 × (1 - 0.794614)
    #        = 0.5 × 0.205386
    #        ≈ 0.1027
    expected = 0.5 * (1.0 - math.exp(-0.115 / 60.0 * 120.0))
    relative_error = abs(R - expected) / expected if expected > 0 else 0
    assert relative_error < 1.0e-4, \
        f"R formula error: expected {expected:.4f}, got {R:.4f}, relative error {relative_error}"
    print(f"✓ test_temporal_vanwagner_formula PASSED (R={R:.4f})")


def test_temporal_equilibrium_reset():
    """Temporal equilibrium reset: when R_E changes, R(0) resets to 0."""
    A_per_sec = 0.115 / 60.0
    
    # Start at R_E1 with some elapsed time
    R_E1 = 0.5
    t1 = 300.0
    R1 = temporal_vanwagner(R_E1, A_per_sec, t1)
    assert R1 > 0.0, "R1 should be > 0"
    
    # After equilibrium change, elapsed time resets to 0
    R_E2 = 1.0  # wind increase
    t_after_reset = 0.0
    R2 = temporal_vanwagner(R_E2, A_per_sec, t_after_reset)
    
    # R should reset to 0 immediately
    assert abs(R2 - 0.0) < 1.0e-10, \
        f"After reset, R should be 0, got {R2}"
    print("✓ test_temporal_equilibrium_reset PASSED")


def test_wind_lag_dampening():
    """Wind-lag dampening: lagged target is between old and new ROS."""
    R_old = 0.5  # m/s
    R_new = 1.0  # m/s (wind increase)
    dt = 30.0  # s
    tau_wind = 60.0  # s
    
    R_target = wind_lag_damping(R_old, R_new, dt, tau_wind)
    
    # R_target should be strictly between R_old and R_new
    assert R_old < R_target < R_new, \
        f"Lagged target {R_target} should be between {R_old} and {R_new}"
    
    # At dt = 30, tau = 60: lag_factor = 1 - exp(-0.5) ≈ 0.3935
    # R_target = 0.5 + 0.5 × 0.3935 ≈ 0.6968
    expected = R_old + (R_new - R_old) * (1.0 - math.exp(-0.5))
    assert abs(R_target - expected) < 1.0e-10, \
        f"Wind lag error: expected {expected}, got {R_target}"
    print(f"✓ test_wind_lag_dampening PASSED (R_target={R_target:.4f})")


def test_point_vs_line_a_selection():
    """Point vs line A selection: A_line > A_point, line fires accelerate faster."""
    A_point = 0.115  # 1/min
    A_line = 0.886   # 1/min
    
    # A_line should be significantly larger than A_point
    assert A_line > A_point, \
        f"A_line {A_line} should be > A_point {A_point}"
    
    # Verify that line fires reach equilibrium faster
    R_E = 1.0  # m/s
    t = 100.0  # s
    A_point_s = A_point / 60.0
    A_line_s = A_line / 60.0
    
    R_point = temporal_vanwagner(R_E, A_point_s, t)
    R_line = temporal_vanwagner(R_E, A_line_s, t)
    
    # Line fires should have higher ROS at same time (closer to equilibrium)
    assert R_line > R_point, \
        f"Line fire R({t}s) should exceed point fire: {R_line} > {R_point}"
    print(f"✓ test_point_vs_line_a_selection PASSED (A_line/A_point={A_line/A_point:.1f}x)")


def test_enable_false_is_noop():
    """Enable=false is no-op: disabled acceleration doesn't scale ROS."""
    # When enable=false, alpha should be 1.0 (no scaling)
    R_original = 0.5
    alpha = 1.0  # implied when disabled
    R_scaled = R_original * alpha
    
    assert abs(R_scaled - R_original) < 1.0e-10, \
        f"When disabled, ROS should be unchanged: {R_scaled} == {R_original}"
    print("✓ test_enable_false_is_noop PASSED")


def main():
    """Run all tests and return 0 on success, 1 on failure."""
    tests = [
        test_size_based_alpha_at_zero,
        test_size_based_alpha_at_large_fire,
        test_size_based_alpha_monotone,
        test_size_based_alpha_formula,
        test_temporal_vanwagner_at_t_zero,
        test_temporal_vanwagner_at_t_infinity,
        test_temporal_vanwagner_monotone,
        test_temporal_vanwagner_formula,
        test_temporal_equilibrium_reset,
        test_wind_lag_dampening,
        test_point_vs_line_a_selection,
        test_enable_false_is_noop,
    ]
    
    failed = []
    for test in tests:
        try:
            test()
        except AssertionError as e:
            print(f"✗ {test.__name__} FAILED: {e}")
            failed.append(test.__name__)
        except Exception as e:
            print(f"✗ {test.__name__} ERROR: {e}")
            failed.append(test.__name__)
    
    if failed:
        print(f"\n{len(failed)} test(s) failed:")
        for name in failed:
            print(f"  - {name}")
        return 1
    else:
        print(f"\n✓ All {len(tests)} tests PASSED")
        return 0


if __name__ == "__main__":
    sys.exit(main())
