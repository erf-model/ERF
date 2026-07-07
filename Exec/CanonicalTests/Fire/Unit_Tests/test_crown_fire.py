#!/usr/bin/env python3
"""Unit tests for crown fire formulas and flame diagnostics."""

import math
import sys

LARGE_INTENSITY = 1.0e12
SIGMA = 5.67e-8
CP_AIR = 1005.0
G = 9.81
PI = math.pi


def van_wagner_critical_intensity(cbh, foliar_moisture, m_c):
    if foliar_moisture <= m_c:
        return LARGE_INTENSITY
    fmc_excess_pct = max((foliar_moisture - m_c) * 100.0, 0.0)
    return 0.010 * max(cbh, 0.1) * (460.0 + 25.9 * fmc_excess_pct)


def cruz_crown_ros(u10_ms, cbd, moisture_10hr):
    if u10_ms <= 0.0:
        return 0.0
    u_kmh = max(u10_ms * 3.6, 1.0)
    mc10_pct = max(moisture_10hr * 100.0, 0.0)
    r_m_min = 11.02 * (u_kmh ** 0.90) * (max(cbd, 0.01) ** 0.19) * math.exp(-0.17 * mc10_pct)
    return max(r_m_min / 60.0, 0.0)


def compute_rothermel_1991_crown_ros(r_surface):
    return max(3.34 * r_surface, 0.0)


def compute_van_wagner_passive_blend(r_surface, r_active, i_b, i_o):
    if i_o < 1.0e-6 or i_b <= 0.0:
        return r_surface
    cf = min(max(i_b / i_o, 0.0) ** (2.0 / 3.0), 1.0)
    return r_surface * (1.0 - cf) + max(r_active, r_surface) * cf


def compute_crown_fraction_burned(r_total, r_surface, r_active):
    r_total = max(r_total, 0.0)
    r_surface = max(r_surface, 0.0)
    r_active = max(r_active, 0.0)
    if r_active <= r_surface + 1.0e-8:
        return 0.0
    return max(0.0, min((r_total - r_surface) / (r_active - r_surface), 1.0))


def compute_flame_length(i_b):
    return 0.0 if i_b <= 0.0 else max(0.0, 0.0775 * (i_b ** 0.46))


def compute_flame_temp_byram(i_b, t_amb):
    return 0.0 if i_b <= 0.0 else max(0.0, t_amb + 800.0 * ((i_b / 1000.0) ** 0.25))


def compute_flame_temp_nelson(i_b):
    if i_b <= 0.0:
        return 0.0
    t_flame = (i_b * 1000.0 / (SIGMA * 0.9)) ** 0.25
    return max(300.0, min(t_flame, 1800.0))


def compute_flame_tilt_angle_deg(i_b, wind_speed, rho_air, t_amb):
    if i_b <= 0.0 or wind_speed <= 0.0 or rho_air <= 0.0 or t_amb <= 0.0:
        return 0.0
    flame_length = compute_flame_length(i_b)
    t_flame = compute_flame_temp_byram(i_b, t_amb)
    q_rad = max(SIGMA * (t_flame ** 4 - t_amb ** 4), 0.0)
    if flame_length <= 0.0 or q_rad <= 0.0:
        return 0.0
    buoyancy = math.sqrt(G * q_rad * flame_length / (rho_air * CP_AIR * t_amb))
    if buoyancy <= 0.0:
        return 90.0
    return min(math.atan(wind_speed / buoyancy) * 180.0 / PI, 90.0)


def crown_initiates(i_b, cbh=2.0, foliar_moisture=0.50, m_c=0.30):
    return i_b >= van_wagner_critical_intensity(cbh, foliar_moisture, m_c)


def test_van_wagner_formula():
    expected = 0.010 * 5.0 * (460.0 + 25.9 * 70.0)
    actual = van_wagner_critical_intensity(5.0, 1.0, 0.30)
    assert abs(actual - expected) < 1.0e-10, f"expected {expected}, got {actual}"
    print("✓ test_van_wagner_formula PASSED")


def test_van_wagner_increases_with_height():
    low = van_wagner_critical_intensity(2.0, 1.0, 0.30)
    high = van_wagner_critical_intensity(8.0, 1.0, 0.30)
    assert high > low, f"higher CBH should increase I_B_crit: {high} <= {low}"
    print("✓ test_van_wagner_increases_with_height PASSED")


def test_van_wagner_impossible_below_Mc():
    actual = van_wagner_critical_intensity(5.0, 0.25, 0.30)
    assert actual >= LARGE_INTENSITY, f"expected impossible threshold, got {actual}"
    print("✓ test_van_wagner_impossible_below_Mc PASSED")


def test_van_wagner_at_Mc_boundary():
    actual = van_wagner_critical_intensity(5.0, 0.30, 0.30)
    assert actual >= LARGE_INTENSITY, f"expected boundary to remain impossible, got {actual}"
    print("✓ test_van_wagner_at_Mc_boundary PASSED")


def test_cruz_formula():
    ros = cruz_crown_ros(30.0 / 3.6, 0.2, 0.10) * 60.0
    assert 15.0 < ros < 45.0, f"expected 15-45 m/min, got {ros}"
    print("✓ test_cruz_formula PASSED")


def test_cruz_increases_with_wind():
    slow = cruz_crown_ros(10.0 / 3.6, 0.2, 0.10)
    fast = cruz_crown_ros(40.0 / 3.6, 0.2, 0.10)
    assert fast > slow, f"higher wind should increase ROS: {fast} <= {slow}"
    print("✓ test_cruz_increases_with_wind PASSED")


def test_cruz_decreases_with_moisture():
    dry = cruz_crown_ros(30.0 / 3.6, 0.2, 0.06)
    wet = cruz_crown_ros(30.0 / 3.6, 0.2, 0.14)
    assert dry > wet, f"higher moisture should reduce ROS: {dry} <= {wet}"
    print("✓ test_cruz_decreases_with_moisture PASSED")


def test_cruz_increases_with_CBD():
    sparse = cruz_crown_ros(30.0 / 3.6, 0.10, 0.10)
    dense = cruz_crown_ros(30.0 / 3.6, 0.30, 0.10)
    assert dense > sparse, f"higher CBD should increase ROS: {dense} <= {sparse}"
    print("✓ test_cruz_increases_with_CBD PASSED")


def test_crown_initiation_logic():
    assert crown_initiates(50.0), "I_B=50 should trigger crown initiation"
    assert not crown_initiates(5.0), "I_B=5 should remain surface fire"
    print("✓ test_crown_initiation_logic PASSED")


def test_cruz_zero_wind():
    ros = cruz_crown_ros(0.0, 0.2, 0.10)
    assert ros == 0.0, f"expected zero ROS at zero wind, got {ros}"
    print("✓ test_cruz_zero_wind PASSED")


def test_rothermel1991_crown_ros():
    ros = compute_rothermel_1991_crown_ros(2.0)
    assert abs(ros - 6.68) < 1.0e-12, f"expected 6.68, got {ros}"
    print("✓ test_rothermel1991_crown_ros PASSED")


def test_passive_blend_at_threshold():
    at_threshold = compute_van_wagner_passive_blend(0.0, 1.0, 10.0, 10.0)
    half_threshold = compute_van_wagner_passive_blend(0.0, 1.0, 5.0, 10.0)
    assert abs(at_threshold - 1.0) < 1.0e-12, f"expected CF=1.0, got {at_threshold}"
    assert 0.60 < half_threshold < 0.66, f"expected CF≈0.63, got {half_threshold}"
    print("✓ test_passive_blend_at_threshold PASSED")


def test_cfb_zero_no_crown():
    cfb = compute_crown_fraction_burned(1.0, 1.0, 3.0)
    assert cfb == 0.0, f"expected 0.0, got {cfb}"
    print("✓ test_cfb_zero_no_crown PASSED")


def test_cfb_one_active_crown():
    cfb = compute_crown_fraction_burned(3.0, 1.0, 3.0)
    assert abs(cfb - 1.0) < 1.0e-12, f"expected 1.0, got {cfb}"
    print("✓ test_cfb_one_active_crown PASSED")


def test_flame_tilt_vertical_no_wind():
    tilt = compute_flame_tilt_angle_deg(500.0, 0.0, 1.2, 300.0)
    assert abs(tilt) < 1.0e-12, f"expected ~0°, got {tilt}"
    print("✓ test_flame_tilt_vertical_no_wind PASSED")


def test_flame_tilt_increases_with_wind():
    low = compute_flame_tilt_angle_deg(500.0, 1.0, 1.2, 300.0)
    high = compute_flame_tilt_angle_deg(500.0, 5.0, 1.2, 300.0)
    assert high > low, f"higher wind should increase tilt: {high} <= {low}"
    print("✓ test_flame_tilt_increases_with_wind PASSED")


def test_byram_radiant_positive():
    t_flame = compute_flame_temp_byram(500.0, 300.0)
    assert 300.0 < t_flame < 1800.0, f"expected 300<T<1800, got {t_flame}"
    print("✓ test_byram_radiant_positive PASSED")


def test_byram_radiant_monotone():
    cool = compute_flame_temp_byram(200.0, 300.0)
    hot = compute_flame_temp_byram(800.0, 300.0)
    assert hot > cool, f"higher I_B should increase T_flame: {hot} <= {cool}"
    print("✓ test_byram_radiant_monotone PASSED")


def test_nelson_emissivity_range():
    t_flame = compute_flame_temp_nelson(200.0)
    assert 300.0 <= t_flame <= 1800.0, f"expected [300, 1800], got {t_flame}"
    print("✓ test_nelson_emissivity_range PASSED")


def main():
    tests = [
        test_van_wagner_formula,
        test_van_wagner_increases_with_height,
        test_van_wagner_impossible_below_Mc,
        test_van_wagner_at_Mc_boundary,
        test_cruz_formula,
        test_cruz_increases_with_wind,
        test_cruz_decreases_with_moisture,
        test_cruz_increases_with_CBD,
        test_crown_initiation_logic,
        test_cruz_zero_wind,
        test_rothermel1991_crown_ros,
        test_passive_blend_at_threshold,
        test_cfb_zero_no_crown,
        test_cfb_one_active_crown,
        test_flame_tilt_vertical_no_wind,
        test_flame_tilt_increases_with_wind,
        test_byram_radiant_positive,
        test_byram_radiant_monotone,
        test_nelson_emissivity_range,
    ]

    failed = []
    for test in tests:
        try:
            test()
        except AssertionError as exc:
            print(f"✗ {test.__name__} FAILED: {exc}")
            failed.append((test.__name__, str(exc)))

    if failed:
        print(f"{len(failed)} test(s) FAILED")
        return 1
    print(f"All {len(tests)} tests PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
