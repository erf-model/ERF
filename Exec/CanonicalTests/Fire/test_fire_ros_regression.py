#!/usr/bin/env python3
"""
Fire ROS regression test.

Checks ROS values for FM1, FM4 against reference values
computed with the correct Rothermel equations. Reference values
are embedded in this file and should be updated whenever the
fire model physics is intentionally changed.

This test does NOT require running ERF — it is a pure-Python
check of the physics formulas.

Run: python3 test_fire_ros_regression.py
"""

import sys
import math


def compute_rothermel_fm1(moisture_1hr=0.08, wind_ftmin=0.0):
    """Compute Rothermel ROS for FM1 (Short Grass)."""
    w0, sigma, delta, Mx = 0.034, 3500.0, 1.0, 0.12
    h_heat, S_T, S_e, rho_p = 8000.0, 0.0555, 0.010, 32.0
    M_f = moisture_1hr

    w_n = w0 * (1 - S_T)
    rho_b = w0 / delta
    beta = rho_b / rho_p
    
    beta_op = 3.348 * sigma**(-0.8189)
    sigma_1p5 = sigma**1.5
    Gamma_max = sigma_1p5 / (495.0 + 0.0594 * sigma_1p5)
    A = 133.0 * sigma**(-0.7913)
    beta_ratio = beta / beta_op
    Gamma_prime = Gamma_max * (beta_ratio**A) * math.exp(A * (1.0 - beta_ratio))

    rm = min(M_f / Mx, 1.0)
    eta_M = max(0.0, 1 - 2.59*rm + 5.11*rm**2 - 3.52*rm**3)
    eta_s = 0.174 * S_e**(-0.19)

    I_R = Gamma_prime * w_n * h_heat * eta_M * eta_s
    I_R = max(I_R, 0.01)

    xi = math.exp((0.792 + 0.681 * math.sqrt(sigma)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)
    eps_h = math.exp(-138.0 / sigma)
    Q_ig = 250.0 + 1116.0 * M_f
    R0_ft_min = (I_R * xi) / (rho_b * eps_h * Q_ig)

    C = 7.47 * math.exp(-0.8711 * sigma**(-0.55))
    B = 0.02526 * sigma**0.54
    E = 0.715 * math.exp(-3.59e-4 * sigma)
    beta_ratio_E = beta_ratio**(-E)

    phi_w_max = 0.9 * I_R
    U_max_ft_min = 0.0
    if C > 0 and B > 0 and beta_ratio_E > 0:
        U_max_ft_min = (phi_w_max / (C * beta_ratio_E))**(1.0/B)

    U_capped = min(wind_ftmin, U_max_ft_min)
    phi_w = C * (U_capped**B) * beta_ratio_E if U_capped > 0 else 0.0

    ROS_ft_min = R0_ft_min * (1.0 + phi_w)
    ROS_ms = ROS_ft_min * 0.00508

    return {
        'R0_ftmin': R0_ft_min, 'R0_ms': R0_ft_min * 0.00508,
        'I_R': I_R, 'C': C, 'B': B, 'E': E,
        'beta': beta, 'beta_op': beta_op, 'beta_ratio_E': beta_ratio_E,
        'U_max_ftmin': U_max_ft_min, 'phi_w': phi_w,
        'ROS_ftmin': ROS_ft_min, 'ROS_ms': ROS_ms
    }


def compute_rothermel_fm4(moisture_1hr=0.08, wind_ftmin=0.0):
    """Compute Rothermel ROS for FM4 (Chaparral)."""
    w0, sigma, delta, Mx = 0.736, 1739.0, 6.0, 0.20
    h_heat, S_T, S_e, rho_p = 8000.0, 0.0555, 0.010, 32.0
    M_f = moisture_1hr

    w_n = w0 * (1 - S_T)
    rho_b = w0 / delta
    beta = rho_b / rho_p
    
    beta_op = 3.348 * sigma**(-0.8189)
    sigma_1p5 = sigma**1.5
    Gamma_max = sigma_1p5 / (495.0 + 0.0594 * sigma_1p5)
    A = 133.0 * sigma**(-0.7913)
    beta_ratio = beta / beta_op
    Gamma_prime = Gamma_max * (beta_ratio**A) * math.exp(A * (1.0 - beta_ratio))

    rm = min(M_f / Mx, 1.0)
    eta_M = max(0.0, 1 - 2.59*rm + 5.11*rm**2 - 3.52*rm**3)
    eta_s = 0.174 * S_e**(-0.19)

    I_R = Gamma_prime * w_n * h_heat * eta_M * eta_s
    I_R = max(I_R, 0.01)

    xi = math.exp((0.792 + 0.681 * math.sqrt(sigma)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)
    eps_h = math.exp(-138.0 / sigma)
    Q_ig = 250.0 + 1116.0 * M_f
    R0_ft_min = (I_R * xi) / (rho_b * eps_h * Q_ig)

    C = 7.47 * math.exp(-0.8711 * sigma**(-0.55))
    B = 0.02526 * sigma**0.54
    E = 0.715 * math.exp(-3.59e-4 * sigma)
    beta_ratio_E = beta_ratio**(-E)

    phi_w_max = 0.9 * I_R
    U_max_ft_min = 0.0
    if C > 0 and B > 0 and beta_ratio_E > 0:
        U_max_ft_min = (phi_w_max / (C * beta_ratio_E))**(1.0/B)

    U_capped = min(wind_ftmin, U_max_ft_min)
    phi_w = C * (U_capped**B) * beta_ratio_E if U_capped > 0 else 0.0

    ROS_ft_min = R0_ft_min * (1.0 + phi_w)
    ROS_ms = ROS_ft_min * 0.00508

    return {
        'R0_ftmin': R0_ft_min, 'R0_ms': R0_ft_min * 0.00508,
        'I_R': I_R, 'C': C, 'B': B, 'E': E,
        'beta': beta, 'beta_op': beta_op, 'beta_ratio_E': beta_ratio_E,
        'U_max_ftmin': U_max_ft_min, 'phi_w': phi_w,
        'ROS_ftmin': ROS_ft_min, 'ROS_ms': ROS_ms
    }


# Reference values (computed with correct Rothermel Eq. 47: C = 7.47*exp(-0.8711*sigma^(-0.55)))
# Units: ROS in m/s, I_R in BTU/ft²/min
# These are expected ranges; actual values should fall within these bounds
REFERENCE_VALUES = {
    # FM1, 8% moisture, no wind
    ('FM1', 0.08, 0.0): {
        'R0_ftmin': (1.5, 5.0),
        'I_R': (100, 1000),
        'C': (7.0, 7.5),
        'U_max_ftmin': (1.0, 100),
    },
    # FM1, 8% moisture, 500 ft/min midflame wind
    ('FM1', 0.08, 500.0): {
        'ROS_ms': (1.0, 50.0),  # Wind significantly increases ROS
        'phi_w': (100, 10000),
    },
    # FM4, 8% moisture, no wind
    ('FM4', 0.08, 0.0): {
        'R0_ftmin': (1.0, 50.0),
        'I_R': (500, 50000),
        'C': (7.0, 7.5),
    },
}


def test_reference_values():
    """Check that computed values match reference ranges."""
    
    print("=" * 70)
    print("Fire ROS Regression Test (Reference Values)")
    print("=" * 70)
    
    failures = []
    test_count = 0
    
    # Test FM1, no wind
    print("\nTesting FM1 at 8% moisture, no wind...")
    r = compute_rothermel_fm1(moisture_1hr=0.08, wind_ftmin=0.0)
    ref = REFERENCE_VALUES[('FM1', 0.08, 0.0)]
    
    for key, (min_val, max_val) in ref.items():
        test_count += 1
        val = r[key]
        if min_val <= val <= max_val:
            print(f"  ✓ {key}: {val:.4f} (expected {min_val:.4f}-{max_val:.4f})")
        else:
            msg = f"{key}: {val:.4f} (expected {min_val:.4f}-{max_val:.4f})"
            print(f"  ✗ {msg}")
            failures.append(msg)
    
    # Test FM1 with wind
    print("\nTesting FM1 at 8% moisture, 500 ft/min wind...")
    r = compute_rothermel_fm1(moisture_1hr=0.08, wind_ftmin=500.0)
    ref = REFERENCE_VALUES[('FM1', 0.08, 500.0)]
    
    for key, (min_val, max_val) in ref.items():
        test_count += 1
        val = r[key]
        if min_val <= val <= max_val:
            print(f"  ✓ {key}: {val:.4f} (expected {min_val:.4f}-{max_val:.4f})")
        else:
            msg = f"{key}: {val:.4f} (expected {min_val:.4f}-{max_val:.4f})"
            print(f"  ✗ {msg}")
            failures.append(msg)
    
    # Test FM4, no wind
    print("\nTesting FM4 at 8% moisture, no wind...")
    r = compute_rothermel_fm4(moisture_1hr=0.08, wind_ftmin=0.0)
    ref = REFERENCE_VALUES[('FM4', 0.08, 0.0)]
    
    for key, (min_val, max_val) in ref.items():
        test_count += 1
        val = r[key]
        if min_val <= val <= max_val:
            print(f"  ✓ {key}: {val:.4f} (expected {min_val:.4f}-{max_val:.4f})")
        else:
            msg = f"{key}: {val:.4f} (expected {min_val:.4f}-{max_val:.4f})"
            print(f"  ✗ {msg}")
            failures.append(msg)
    
    print("\n" + "=" * 70)
    if not failures:
        print(f"FINAL RESULT: All {test_count} reference value checks PASSED")
        return 0
    else:
        print(f"FINAL RESULT: {len(failures)} of {test_count} checks FAILED")
        for failure in failures:
            print(f"  ✗ {failure}")
        return 1


if __name__ == "__main__":
    sys.exit(test_reference_values())
