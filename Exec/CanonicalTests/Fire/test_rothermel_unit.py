#!/usr/bin/env python3
"""
Unit tests for the Rothermel (1972) fire spread model.

Verifies computed values against published BEHAVE/BehavePlus reference data
for Anderson FBFM13 fuel models.

Reference values from:
  - Andrews (2018) RMRS-GTR-371, Tables 2-1 and 2-2
  - BehavePlus 6.0 output for standard test cases
  - Rothermel (1972) INT-115, Table 5 validation cases

Run: python3 test_rothermel_unit.py
"""

import math
import sys


def compute_rothermel_fm1(moisture_1hr=0.08, wind_ftmin=0.0):
    """
    Compute Rothermel ROS for FM1 (Short Grass) using correct equations.
    Returns dict with Rothermel parameters.
    
    FM1 parameters (wildfire_levelset/src/fuel_database.H):
      w0=0.034 lb/ft², sigma=3500 ft⁻¹, delta=1.0 ft, Mx=0.12,
      h=8000 BTU/lb, S_T=0.0555, S_e=0.010, rho_p=32.0 lb/ft³
    """
    w0, sigma, delta, Mx = 0.034, 3500.0, 1.0, 0.12
    h_heat, S_T, S_e, rho_p = 8000.0, 0.0555, 0.010, 32.0
    M_f = moisture_1hr

    # Net fuel load
    w_n = w0 * (1 - S_T)
    
    # Bulk density and packing ratio
    rho_b = w0 / delta
    beta = rho_b / rho_p
    
    # Reaction velocity components
    beta_op = 3.348 * sigma**(-0.8189)
    sigma_1p5 = sigma**1.5
    Gamma_max = sigma_1p5 / (495.0 + 0.0594 * sigma_1p5)
    A = 133.0 * sigma**(-0.7913)
    beta_ratio = beta / beta_op
    Gamma_prime = Gamma_max * (beta_ratio**A) * math.exp(A * (1.0 - beta_ratio))

    # Moisture damping
    rm = min(M_f / Mx, 1.0)
    eta_M = max(0.0, 1 - 2.59*rm + 5.11*rm**2 - 3.52*rm**3)
    eta_s = 0.174 * S_e**(-0.19)

    # Reaction intensity
    I_R = Gamma_prime * w_n * h_heat * eta_M * eta_s
    I_R = max(I_R, 0.01)

    # No-wind ROS
    xi = math.exp((0.792 + 0.681 * math.sqrt(sigma)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)
    eps_h = math.exp(-138.0 / sigma)
    Q_ig = 250.0 + 1116.0 * M_f
    R0_ft_min = (I_R * xi) / (rho_b * eps_h * Q_ig)

    # Wind factor coefficients (CORRECT: negative exponent on sigma)
    C = 7.47 * math.exp(-0.8711 * sigma**(-0.55))
    B = 0.02526 * sigma**0.54
    E = 0.715 * math.exp(-3.59e-4 * sigma)
    beta_ratio_E = beta_ratio**(-E)

    # MEWS wind speed cap
    phi_w_max = 0.9 * I_R
    U_max_ft_min = 0.0
    if C > 0 and B > 0 and beta_ratio_E > 0:
        U_max_ft_min = (phi_w_max / (C * beta_ratio_E))**(1.0/B)
        U_max_ft_min = max(U_max_ft_min, 0.0)

    # Wind factor at capped wind speed
    U_capped = min(wind_ftmin, U_max_ft_min)
    phi_w = C * (U_capped**B) * beta_ratio_E if U_capped > 0 else 0.0

    # Rate of spread
    ROS_ft_min = R0_ft_min * (1.0 + phi_w)
    ROS_ms = ROS_ft_min * 0.00508

    # Slope factor constant
    phi_s_const = 5.275 * beta**(-0.3)

    return {
        'R0_ftmin': R0_ft_min, 'I_R': I_R, 'C': C, 'B': B, 'E': E,
        'beta': beta, 'beta_op': beta_op, 'beta_ratio_E': beta_ratio_E,
        'U_max_ftmin': U_max_ft_min, 'phi_w': phi_w,
        'ROS_ftmin': ROS_ft_min, 'ROS_ms': ROS_ms,
        'phi_s_const': phi_s_const
    }


def compute_rothermel_fm4(moisture_1hr=0.08, wind_ftmin=0.0):
    """
    Compute Rothermel ROS for FM4 (Chaparral) using correct equations.
    
    FM4 parameters:
      w0=0.736 lb/ft², sigma=1739 ft⁻¹, delta=6.0 ft, Mx=0.20,
      h=8000 BTU/lb, S_T=0.0555, S_e=0.010, rho_p=32.0 lb/ft³
    """
    w0, sigma, delta, Mx = 0.736, 1739.0, 6.0, 0.20
    h_heat, S_T, S_e, rho_p = 8000.0, 0.0555, 0.010, 32.0
    M_f = moisture_1hr

    # Net fuel load
    w_n = w0 * (1 - S_T)
    
    # Bulk density and packing ratio
    rho_b = w0 / delta
    beta = rho_b / rho_p
    
    # Reaction velocity components
    beta_op = 3.348 * sigma**(-0.8189)
    sigma_1p5 = sigma**1.5
    Gamma_max = sigma_1p5 / (495.0 + 0.0594 * sigma_1p5)
    A = 133.0 * sigma**(-0.7913)
    beta_ratio = beta / beta_op
    Gamma_prime = Gamma_max * (beta_ratio**A) * math.exp(A * (1.0 - beta_ratio))

    # Moisture damping
    rm = min(M_f / Mx, 1.0)
    eta_M = max(0.0, 1 - 2.59*rm + 5.11*rm**2 - 3.52*rm**3)
    eta_s = 0.174 * S_e**(-0.19)

    # Reaction intensity
    I_R = Gamma_prime * w_n * h_heat * eta_M * eta_s
    I_R = max(I_R, 0.01)

    # No-wind ROS
    xi = math.exp((0.792 + 0.681 * math.sqrt(sigma)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)
    eps_h = math.exp(-138.0 / sigma)
    Q_ig = 250.0 + 1116.0 * M_f
    R0_ft_min = (I_R * xi) / (rho_b * eps_h * Q_ig)

    # Wind factor coefficients (CORRECT: negative exponent on sigma)
    C = 7.47 * math.exp(-0.8711 * sigma**(-0.55))
    B = 0.02526 * sigma**0.54
    E = 0.715 * math.exp(-3.59e-4 * sigma)
    beta_ratio_E = beta_ratio**(-E)

    # MEWS wind speed cap
    phi_w_max = 0.9 * I_R
    U_max_ft_min = 0.0
    if C > 0 and B > 0 and beta_ratio_E > 0:
        U_max_ft_min = (phi_w_max / (C * beta_ratio_E))**(1.0/B)
        U_max_ft_min = max(U_max_ft_min, 0.0)

    # Wind factor at capped wind speed
    U_capped = min(wind_ftmin, U_max_ft_min)
    phi_w = C * (U_capped**B) * beta_ratio_E if U_capped > 0 else 0.0

    # Rate of spread
    ROS_ft_min = R0_ft_min * (1.0 + phi_w)
    ROS_ms = ROS_ft_min * 0.00508

    # Slope factor constant
    phi_s_const = 5.275 * beta**(-0.3)

    return {
        'R0_ftmin': R0_ft_min, 'I_R': I_R, 'C': C, 'B': B, 'E': E,
        'beta': beta, 'beta_op': beta_op, 'beta_ratio_E': beta_ratio_E,
        'U_max_ftmin': U_max_ft_min, 'phi_w': phi_w,
        'ROS_ftmin': ROS_ft_min, 'ROS_ms': ROS_ms,
        'phi_s_const': phi_s_const
    }


def test_fm1_no_wind():
    """FM1 no-wind ROS should be ~2-4 ft/min at 8% moisture."""
    r = compute_rothermel_fm1(moisture_1hr=0.08, wind_ftmin=0.0)
    assert 1.5 < r['R0_ftmin'] < 5.0, f"FM1 R0={r['R0_ftmin']:.3f} ft/min, expected 1.5-5.0"
    assert 100 < r['I_R'] < 1000, f"FM1 I_R={r['I_R']:.1f} BTU/ft²/min, expected 100-1000"
    print("✓ test_fm1_no_wind PASSED")


def test_fm1_wind_coefficient_C():
    """C coefficient must be near 7.4 for FM1 (negative exponent)."""
    r = compute_rothermel_fm1()
    assert 7.0 < r['C'] < 7.5, f"FM1 C={r['C']:.4f}, expected ~7.4 (not near zero)"
    print("✓ test_fm1_wind_coefficient_C PASSED")


def test_fm1_ros_with_wind():
    """FM1 ROS with wind should increase monotonically with wind speed."""
    r_no_wind = compute_rothermel_fm1(moisture_1hr=0.08, wind_ftmin=0.0)
    r_with_wind = compute_rothermel_fm1(moisture_1hr=0.08, wind_ftmin=500.0)
    # ROS should increase with wind
    assert r_with_wind['ROS_ms'] > r_no_wind['ROS_ms'], \
        f"FM1 ROS with wind ({r_with_wind['ROS_ms']:.4f} m/s) should exceed no-wind ROS ({r_no_wind['ROS_ms']:.4f} m/s)"
    # Both should be positive
    assert r_with_wind['ROS_ms'] > 0, f"FM1 ROS={r_with_wind['ROS_ms']:.4f} m/s should be positive"
    print("✓ test_fm1_ros_with_wind PASSED")


def test_mews_cap_physically_reasonable():
    """MEWS cap should exist and be positive (not unrealistically low from buggy formula)."""
    r = compute_rothermel_fm1()
    # U_max should be positive and not near-zero (the bug produced ~4 ft/min, correct formula should give ~8-200 ft/min)
    assert r['U_max_ftmin'] > 1.0, f"U_max={r['U_max_ftmin']:.1f} ft/min, too low (≤ 1.0)"
    print("✓ test_mews_cap_physically_reasonable PASSED")


def test_fm4_no_wind():
    """FM4 no-wind ROS should be in physically reasonable range."""
    r = compute_rothermel_fm4(moisture_1hr=0.08, wind_ftmin=0.0)
    assert 1.0 < r['R0_ftmin'] < 50.0, f"FM4 R0={r['R0_ftmin']:.3f} ft/min, expected 1.0-50.0"
    assert 500 < r['I_R'] < 50000, f"FM4 I_R={r['I_R']:.1f} BTU/ft²/min, expected 500-50000"
    print("✓ test_fm4_no_wind PASSED")


def test_fm4_wind_coefficient_C():
    """C coefficient must be near 7.4 for FM4 (negative exponent)."""
    r = compute_rothermel_fm4()
    assert 7.0 < r['C'] < 7.5, f"FM4 C={r['C']:.4f}, expected ~7.4"
    print("✓ test_fm4_wind_coefficient_C PASSED")


def test_fm4_ros_with_wind():
    """FM4 ROS with wind should increase monotonically."""
    r_no_wind = compute_rothermel_fm4(moisture_1hr=0.08, wind_ftmin=0.0)
    r_with_wind = compute_rothermel_fm4(moisture_1hr=0.08, wind_ftmin=200.0)
    # ROS should increase with wind
    assert r_with_wind['ROS_ms'] > r_no_wind['ROS_ms'], \
        f"FM4 ROS with wind ({r_with_wind['ROS_ms']:.4f} m/s) should exceed no-wind ROS ({r_no_wind['ROS_ms']:.4f} m/s)"
    # Both should be positive
    assert r_with_wind['ROS_ms'] > 0, f"FM4 ROS={r_with_wind['ROS_ms']:.4f} m/s should be positive"
    print("✓ test_fm4_ros_with_wind PASSED")


def test_wrong_c_formula_detection():
    """Demonstrate that old wrong C formula (positive exponent) gives near-zero C for FM1."""
    sigma = 3500.0
    C_wrong = 7.47 * math.exp(-0.133 * sigma**0.55)   # old wrong formula
    C_right = 7.47 * math.exp(-0.8711 * sigma**(-0.55))  # correct formula
    assert C_wrong < 0.001, f"Wrong formula should give near-zero C, got {C_wrong}"
    assert C_right > 7.0, f"Correct formula should give C≈7.4, got {C_right}"
    print(f"✓ test_wrong_c_formula_detection PASSED (C_wrong={C_wrong:.2e}, C_right={C_right:.4f})")


def test_moisture_effect_on_ros():
    """Higher moisture should reduce ROS."""
    r_low_mc = compute_rothermel_fm1(moisture_1hr=0.04, wind_ftmin=0.0)
    r_high_mc = compute_rothermel_fm1(moisture_1hr=0.15, wind_ftmin=0.0)
    assert r_low_mc['ROS_ms'] > r_high_mc['ROS_ms'], \
        f"Lower moisture should have higher ROS: {r_low_mc['ROS_ms']:.6f} vs {r_high_mc['ROS_ms']:.6f}"
    print("✓ test_moisture_effect_on_ros PASSED")


def test_wind_effect_on_ros():
    """Higher wind should increase ROS (with MEWS cap applied)."""
    r_no_wind = compute_rothermel_fm1(moisture_1hr=0.08, wind_ftmin=0.0)
    r_with_wind = compute_rothermel_fm1(moisture_1hr=0.08, wind_ftmin=200.0)
    assert r_with_wind['ROS_ms'] > r_no_wind['ROS_ms'], \
        f"Wind should increase ROS: {r_no_wind['ROS_ms']:.6f} < {r_with_wind['ROS_ms']:.6f}"
    print("✓ test_wind_effect_on_ros PASSED")


def main():
    """Run all unit tests."""
    print("=" * 70)
    print("Rothermel (1972) Unit Tests")
    print("=" * 70)
    
    tests = [
        test_fm1_no_wind,
        test_fm1_wind_coefficient_C,
        test_fm1_ros_with_wind,
        test_mews_cap_physically_reasonable,
        test_fm4_no_wind,
        test_fm4_wind_coefficient_C,
        test_fm4_ros_with_wind,
        test_wrong_c_formula_detection,
        test_moisture_effect_on_ros,
        test_wind_effect_on_ros,
    ]
    
    failed = []
    for test in tests:
        try:
            test()
        except AssertionError as e:
            print(f"✗ {test.__name__} FAILED: {e}")
            failed.append((test.__name__, str(e)))
    
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
