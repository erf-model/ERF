#!/usr/bin/env python3
"""
Unit tests for Phase 13 rate-of-spread models: MacArthur, Balbi, Cheney-Gould, 
and per-fuel wind height functions.

Verifies pure-Python implementations of ROS model formulae against analytical
reference values extracted from ERF source code.

Tests are organized by model:
  - MacArthur (1966) Australian formula (4 tests)
  - Balbi (2009) physical model (5 tests)
  - Cheney-Gould (1998) grassland model (4 tests)
  - Per-fuel wind height tables (5 tests)

Reference implementations extracted from:
  - ERF_BalbiModel.H: macarthur_ros(), compute_balbi_angle(), ROS formula
  - ERF_CheneyGouldModel.H: cheney_gould_ros()
  - ERF_FuelWindHeight.H: build_fcwh_table(), build_fcz0_table()

Run: python3 test_ros_models.py
"""

import math
import sys


# ============================================================================
# MacArthur (1966) Australian Formula Tests
# ============================================================================

def macarthur_ros(wind_speed_ms):
    """
    MacArthur (1966) Mark 5 Forest Fire Danger Meter ROS formula.
    
    Formula from ERF_BalbiModel.H lines 159-165:
      R [m/s] = backing * exp(0.8424 * max(U, 0))
      where backing = 0.18 m/s
    
    Args:
        wind_speed_ms: Wind speed in fire spread direction [m/s]
    
    Returns:
        ROS [m/s]
    """
    backing = 0.18
    wind_contribution = backing * (math.exp(0.8424 * max(wind_speed_ms, 0.0)) - 1.0)
    return backing + wind_contribution


def test_macarthur_no_wind_backing():
    """Test 1: MacArthur no-wind backing rate equals 0.18 m/s."""
    ros = macarthur_ros(0.0)
    expected = 0.18
    passed = abs(ros - expected) < 1e-6
    status = "✓" if passed else "✗"
    print(f"{status} Test 1: MacArthur no-wind backing rate")
    if not passed:
        print(f"    Expected: {expected}, Got: {ros}")
    return passed


def test_macarthur_formula_structure():
    """Test 2: MacArthur formula structure at U=5 m/s."""
    U = 5.0
    ros = macarthur_ros(U)
    # Verify formula: R = 0.18 + 0.18 * (exp(0.8424*5) - 1)
    expected = 0.18 + 0.18 * (math.exp(0.8424 * 5) - 1.0)
    passed = abs(ros - expected) < 1e-6 and ros > 0.18
    status = "✓" if passed else "✗"
    print(f"{status} Test 2: MacArthur formula structure at U=5 m/s")
    if not passed:
        print(f"    Expected: {expected}, Got: {ros}, ros > 0.18: {ros > 0.18}")
    return passed


def test_macarthur_monotone_with_wind():
    """Test 3: MacArthur ROS increases monotonically with wind."""
    ros_0 = macarthur_ros(0.0)
    ros_2 = macarthur_ros(2.0)
    ros_5 = macarthur_ros(5.0)
    passed = (ros_0 < ros_2) and (ros_2 < ros_5)
    status = "✓" if passed else "✗"
    print(f"{status} Test 3: MacArthur ROS increases monotonically with wind")
    if not passed:
        print(f"    R(0)={ros_0}, R(2)={ros_2}, R(5)={ros_5}")
    return passed


def test_macarthur_negative_wind_clamp():
    """Test 4: MacArthur clamps negative wind to backing rate."""
    # Opposing wind (U = -3 m/s) clamps to 0, so R = backing = 0.18
    ros = macarthur_ros(-3.0)
    expected = 0.18
    passed = abs(ros - expected) < 1e-6
    status = "✓" if passed else "✗"
    print(f"{status} Test 4: MacArthur negative wind clamp to backing rate")
    if not passed:
        print(f"    Expected: {expected}, Got: {ros}")
    return passed


# ============================================================================
# Balbi (2009) Physical Model Tests
# ============================================================================

def compute_balbi_angle(U, theta, v_b):
    """
    Compute wind and slope angle for Balbi model.
    
    Formula from ERF_BalbiModel.H lines 127-137:
      tan α = U/v_b + tan θ
      where v_b < 1e-6 is clamped to 1e-6 to avoid division by zero
    
    Args:
        U: Wind speed in spread direction [m/s]
        theta: Terrain slope [radians]
        v_b: Buoyancy velocity [m/s]
    
    Returns:
        Angle α [radians]
    """
    small_v_b = 1.0e-6
    v_b = max(v_b, small_v_b)
    tan_alpha = U / v_b + math.tan(theta)
    return math.atan(tan_alpha)


def balbi_ros(U, theta, A_coeff, v_b):
    """
    Balbi (2009) ROS formula.
    
    Formula from ERF_BalbiModel.H lines 233-238:
      R = A_coeff × (1 + sin α − cos α)
      where α = compute_balbi_angle(U, theta, v_b)
    
    Args:
        U: Wind speed in spread direction [m/s]
        theta: Terrain slope [radians]
        A_coeff: Balbi coefficient A [m/s]
        v_b: Buoyancy velocity [m/s]
    
    Returns:
        ROS [m/s], clamped to ≥ 0
    """
    alpha = compute_balbi_angle(U, theta, v_b)
    sin_alpha = math.sin(alpha)
    cos_alpha = math.cos(alpha)
    ros = A_coeff * (1.0 + sin_alpha - cos_alpha)
    return max(ros, 0.0)


def test_balbi_zero_wind_slope():
    """Test 5: Balbi returns zero ROS at zero wind and zero slope."""
    # At U=0, θ=0: tan α = 0/v_b + 0 = 0, so α = 0
    # R = A × (1 + 0 - 1) = 0
    A_coeff = 0.5
    v_b = 2.0
    ros = balbi_ros(U=0.0, theta=0.0, A_coeff=A_coeff, v_b=v_b)
    expected = 0.0
    passed = abs(ros - expected) < 1e-6
    status = "✓" if passed else "✗"
    print(f"{status} Test 5: Balbi zero ROS at zero wind and zero slope")
    if not passed:
        print(f"    Expected: {expected}, Got: {ros}")
    return passed


def test_balbi_increases_with_wind():
    """Test 6: Balbi ROS increases with wind for fixed A_coeff and v_b."""
    A_coeff = 0.5
    v_b = 2.0
    ros_0 = balbi_ros(U=0.0, theta=0.0, A_coeff=A_coeff, v_b=v_b)
    ros_2 = balbi_ros(U=2.0, theta=0.0, A_coeff=A_coeff, v_b=v_b)
    ros_4 = balbi_ros(U=4.0, theta=0.0, A_coeff=A_coeff, v_b=v_b)
    passed = (ros_0 < ros_2) and (ros_2 < ros_4)
    status = "✓" if passed else "✗"
    print(f"{status} Test 6: Balbi ROS increases with wind")
    if not passed:
        print(f"    R(U=0)={ros_0}, R(U=2)={ros_2}, R(U=4)={ros_4}")
    return passed


def test_balbi_angle_formula():
    """Test 7: Balbi angle formula verification."""
    # For U=3 m/s, v_b=2 m/s, slope=0:
    #   tan α = 3/2 = 1.5
    #   α = atan(1.5)
    U = 3.0
    v_b = 2.0
    theta = 0.0
    alpha = compute_balbi_angle(U, theta, v_b)
    expected_alpha = math.atan(1.5)
    
    # Now compute ROS and verify
    A_coeff = 0.5
    ros = balbi_ros(U, theta, A_coeff, v_b)
    sin_alpha = math.sin(expected_alpha)
    cos_alpha = math.cos(expected_alpha)
    expected_ros = A_coeff * (1.0 + sin_alpha - cos_alpha)
    
    passed = abs(alpha - expected_alpha) < 1e-6 and abs(ros - expected_ros) < 1e-6
    status = "✓" if passed else "✗"
    print(f"{status} Test 7: Balbi angle formula (U=3, v_b=2, θ=0)")
    if not passed:
        print(f"    α: expected {expected_alpha}, got {alpha}")
        print(f"    R: expected {expected_ros}, got {ros}")
    return passed


def test_balbi_small_v_b_floor():
    """Test 8: Balbi small v_b floor prevents division by zero."""
    # When v_b < 1e-6, it is clamped to 1e-6
    A_coeff = 0.5
    v_b_tiny = 1.0e-8  # Much smaller than floor
    U = 1.0
    theta = 0.0
    # Should not raise or produce NaN
    ros = balbi_ros(U, theta, A_coeff, v_b_tiny)
    passed = math.isfinite(ros) and ros >= 0.0
    status = "✓" if passed else "✗"
    print(f"{status} Test 8: Balbi small v_b floor prevents division by zero")
    if not passed:
        print(f"    ROS={ros}, is_finite={math.isfinite(ros)}")
    return passed


def test_balbi_positive_ros_clamp():
    """Test 9: Balbi ROS is always non-negative."""
    # Opposing slope with large wind should still produce positive ROS
    A_coeff = 0.5
    v_b = 1.0
    U = 5.0
    theta = -0.2  # Downslope (opposing)
    ros = balbi_ros(U, theta, A_coeff, v_b)
    passed = ros >= 0.0
    status = "✓" if passed else "✗"
    print(f"{status} Test 9: Balbi ROS clamped to non-negative")
    if not passed:
        print(f"    ROS={ros}, ros >= 0: {ros >= 0}")
    return passed


# ============================================================================
# Cheney-Gould (1998) Grassland Model Tests
# ============================================================================

def cheney_gould_ros(u_wind, ros_backing, moisture, curing):
    """
    Cheney-Gould (1998) ROS formula.
    
    Formula from ERF_CheneyGouldModel.H lines 104-114:
      wind_factor = 0.15 * U * (curing + 0.2)
      moisture_factor = 20 / (moisture + 1)
      ROS_forward = ros_backing * (1 + wind_factor) * moisture_factor
    
    Args:
        u_wind: Wind speed in forward direction [m/s]
        ros_backing: Base (no-wind) ROS [m/s]
        moisture: Dead fine fuel moisture [%]
        curing: Curing degree [0-1]
    
    Returns:
        ROS [m/s], clamped to ≥ 0
    """
    wind_factor = 0.15 * u_wind * (curing + 0.2)
    moisture_scale = 20.0
    moisture_factor = moisture_scale / (moisture + 1.0)
    ros_forward = ros_backing * (1.0 + wind_factor) * moisture_factor
    return max(ros_forward, 0.0)


def test_cheney_gould_zero_wind():
    """Test 10: Cheney-Gould zero wind formula."""
    # At u_wind=0:
    # wind_factor = 0.15 * 0 * (curing + 0.2) = 0
    # ROS = ros_backing * (1 + 0) * moisture_factor
    # With ros_backing=0.05, moisture=10, curing=1.0:
    # ROS = 0.05 * 1 * (20/(10+1)) = 0.05 * (20/11)
    ros_backing = 0.05
    moisture = 10.0
    curing = 1.0
    ros = cheney_gould_ros(u_wind=0.0, ros_backing=ros_backing, 
                           moisture=moisture, curing=curing)
    expected = ros_backing * 1.0 * (20.0 / (moisture + 1.0))
    passed = abs(ros - expected) < 1e-6
    status = "✓" if passed else "✗"
    print(f"{status} Test 10: Cheney-Gould zero wind formula")
    if not passed:
        print(f"    Expected: {expected}, Got: {ros}")
    return passed


def test_cheney_gould_increases_with_wind():
    """Test 11: Cheney-Gould ROS increases with wind."""
    ros_backing = 0.05
    moisture = 10.0
    curing = 1.0
    ros_0 = cheney_gould_ros(u_wind=0.0, ros_backing=ros_backing,
                             moisture=moisture, curing=curing)
    ros_3 = cheney_gould_ros(u_wind=3.0, ros_backing=ros_backing,
                             moisture=moisture, curing=curing)
    ros_5 = cheney_gould_ros(u_wind=5.0, ros_backing=ros_backing,
                             moisture=moisture, curing=curing)
    passed = (ros_0 < ros_3) and (ros_3 < ros_5)
    status = "✓" if passed else "✗"
    print(f"{status} Test 11: Cheney-Gould ROS increases with wind")
    if not passed:
        print(f"    R(U=0)={ros_0}, R(U=3)={ros_3}, R(U=5)={ros_5}")
    return passed


def test_cheney_gould_decreases_with_moisture():
    """Test 12: Cheney-Gould ROS decreases with moisture."""
    ros_backing = 0.05
    curing = 1.0
    u_wind = 3.0
    ros_5pct = cheney_gould_ros(u_wind=u_wind, ros_backing=ros_backing,
                                moisture=5.0, curing=curing)
    ros_10pct = cheney_gould_ros(u_wind=u_wind, ros_backing=ros_backing,
                                 moisture=10.0, curing=curing)
    ros_20pct = cheney_gould_ros(u_wind=u_wind, ros_backing=ros_backing,
                                 moisture=20.0, curing=curing)
    passed = (ros_5pct > ros_10pct) and (ros_10pct > ros_20pct)
    status = "✓" if passed else "✗"
    print(f"{status} Test 12: Cheney-Gould ROS decreases with moisture")
    if not passed:
        print(f"    R(5%)={ros_5pct}, R(10%)={ros_10pct}, R(20%)={ros_20pct}")
    return passed


def test_cheney_gould_increases_with_curing():
    """Test 13: Cheney-Gould ROS increases with curing."""
    ros_backing = 0.05
    u_wind = 3.0
    moisture = 10.0
    ros_0_cure = cheney_gould_ros(u_wind=u_wind, ros_backing=ros_backing,
                                  moisture=moisture, curing=0.0)
    ros_0_5_cure = cheney_gould_ros(u_wind=u_wind, ros_backing=ros_backing,
                                    moisture=moisture, curing=0.5)
    ros_1_0_cure = cheney_gould_ros(u_wind=u_wind, ros_backing=ros_backing,
                                    moisture=moisture, curing=1.0)
    passed = (ros_0_cure < ros_0_5_cure) and (ros_0_5_cure < ros_1_0_cure)
    status = "✓" if passed else "✗"
    print(f"{status} Test 13: Cheney-Gould ROS increases with curing")
    if not passed:
        print(f"    R(cure=0)={ros_0_cure}, R(cure=0.5)={ros_0_5_cure}, R(cure=1)={ros_1_0_cure}")
    return passed


# ============================================================================
# Per-Fuel Wind Height Tests
# ============================================================================

def build_fcwh_table(global_z_ref, use_per_fuel=False):
    """
    Build per-fuel wind height (fcwh) table indexed 0..13.
    
    Implementation from ERF_FuelWindHeight.H lines 44-63:
    - When use_per_fuel=False: all entries 1-13 equal global_z_ref
    - When use_per_fuel=True: all entries 1-13 equal 6.096 (WRF-SFIRE default)
    
    Args:
        global_z_ref: Global fallback wind reference height [m]
        use_per_fuel: When True, use WRF-SFIRE defaults; when False, use global_z_ref
    
    Returns:
        List of size 14; index 0 unused, 1-13 valid
    """
    fcwh = [0.0] * 14
    if use_per_fuel:
        for i in range(1, 14):
            fcwh[i] = 6.096
    else:
        for i in range(1, 14):
            fcwh[i] = global_z_ref
    return fcwh


def build_fcz0_table():
    """
    Build per-fuel roughness length (fcz0) table indexed 0..13.
    
    Implementation from ERF_FuelWindHeight.H lines 73-93.
    WRF-SFIRE data statement values [m].
    
    Returns:
        List of size 14; index 0 unused, 1-13 valid
    """
    fcz0 = [0.0] * 14
    fcz0[1]  = 0.0396   # FM1
    fcz0[2]  = 0.0396   # FM2
    fcz0[3]  = 0.100    # FM3
    fcz0[4]  = 0.2378   # FM4
    fcz0[5]  = 0.0793   # FM5
    fcz0[6]  = 0.0991   # FM6
    fcz0[7]  = 0.0991   # FM7
    fcz0[8]  = 0.0079   # FM8
    fcz0[9]  = 0.0079   # FM9
    fcz0[10] = 0.0396   # FM10
    fcz0[11] = 0.0396   # FM11
    fcz0[12] = 0.0911   # FM12
    fcz0[13] = 0.1188   # FM13
    return fcz0


def test_fcwh_uniform_mode():
    """Test 14: fcwh uniform mode returns global_z_ref for all fuels."""
    global_z_ref = 6.1
    fcwh = build_fcwh_table(global_z_ref, use_per_fuel=False)
    passed = (len(fcwh) == 14 and 
              all(fcwh[i] == global_z_ref for i in range(1, 14)))
    status = "✓" if passed else "✗"
    print(f"{status} Test 14: fcwh uniform mode (all fuels = {global_z_ref})")
    if not passed:
        print(f"    Length: {len(fcwh)}, entries 1-13 all {global_z_ref}: "
              f"{all(fcwh[i] == global_z_ref for i in range(1, 14))}")
    return passed


def test_fcwh_per_fuel_mode():
    """Test 15: fcwh per-fuel mode returns 6.096 for all fuels."""
    global_z_ref = 6.1
    fcwh = build_fcwh_table(global_z_ref, use_per_fuel=True)
    expected = 6.096
    passed = (len(fcwh) == 14 and 
              all(abs(fcwh[i] - expected) < 1e-6 for i in range(1, 14)))
    status = "✓" if passed else "✗"
    print(f"{status} Test 15: fcwh per-fuel mode (all fuels = 6.096 m)")
    if not passed:
        print(f"    Length: {len(fcwh)}")
        for i in range(1, 14):
            if abs(fcwh[i] - expected) >= 1e-6:
                print(f"    fcwh[{i}] = {fcwh[i]}, expected {expected}")
    return passed


def test_fcz0_fm4_chaparral():
    """Test 16: fcz0 FM4 value equals 0.2378 (chaparral, highest roughness)."""
    fcz0 = build_fcz0_table()
    expected = 0.2378
    passed = abs(fcz0[4] - expected) < 1e-6
    status = "✓" if passed else "✗"
    print(f"{status} Test 16: fcz0 FM4 = 0.2378 m (chaparral)")
    if not passed:
        print(f"    Expected: {expected}, Got: {fcz0[4]}")
    return passed


def test_fcz0_fm1_fm2_equal():
    """Test 17: fcz0 FM1 and FM2 both equal 0.0396."""
    fcz0 = build_fcz0_table()
    expected = 0.0396
    passed = (abs(fcz0[1] - expected) < 1e-6 and 
              abs(fcz0[2] - expected) < 1e-6 and 
              fcz0[1] == fcz0[2])
    status = "✓" if passed else "✗"
    print(f"{status} Test 17: fcz0 FM1 and FM2 both equal 0.0396 m")
    if not passed:
        print(f"    fcz0[1] = {fcz0[1]}, fcz0[2] = {fcz0[2]}, expected {expected}")
    return passed


def test_fcz0_table_size():
    """Test 18: fcz0 table has size 14 (indices 0-13)."""
    fcz0 = build_fcz0_table()
    passed = len(fcz0) == 14
    status = "✓" if passed else "✗"
    print(f"{status} Test 18: fcz0 table size = 14")
    if not passed:
        print(f"    Expected length 14, got {len(fcz0)}")
    return passed


# ============================================================================
# Test Runner
# ============================================================================

def main():
    """Run all 18 tests and return exit code (0 = all pass, 1 = any fail)."""
    print("=" * 70)
    print("Phase 13 ROS Model Unit Tests")
    print("=" * 70)
    print()
    
    # MacArthur tests (4)
    print("MacArthur (1966) Australian Formula Tests")
    print("-" * 70)
    results = []
    results.append(test_macarthur_no_wind_backing())
    results.append(test_macarthur_formula_structure())
    results.append(test_macarthur_monotone_with_wind())
    results.append(test_macarthur_negative_wind_clamp())
    print()
    
    # Balbi tests (5)
    print("Balbi (2009) Physical Model Tests")
    print("-" * 70)
    results.append(test_balbi_zero_wind_slope())
    results.append(test_balbi_increases_with_wind())
    results.append(test_balbi_angle_formula())
    results.append(test_balbi_small_v_b_floor())
    results.append(test_balbi_positive_ros_clamp())
    print()
    
    # Cheney-Gould tests (4)
    print("Cheney-Gould (1998) Grassland Model Tests")
    print("-" * 70)
    results.append(test_cheney_gould_zero_wind())
    results.append(test_cheney_gould_increases_with_wind())
    results.append(test_cheney_gould_decreases_with_moisture())
    results.append(test_cheney_gould_increases_with_curing())
    print()
    
    # Per-fuel wind height tests (5)
    print("Per-Fuel Wind Height Tests")
    print("-" * 70)
    results.append(test_fcwh_uniform_mode())
    results.append(test_fcwh_per_fuel_mode())
    results.append(test_fcz0_fm4_chaparral())
    results.append(test_fcz0_fm1_fm2_equal())
    results.append(test_fcz0_table_size())
    print()
    
    # Summary
    total = len(results)
    passed = sum(results)
    failed = total - passed
    print("=" * 70)
    print(f"Results: {passed}/{total} passed, {failed}/{total} failed")
    print("=" * 70)
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
