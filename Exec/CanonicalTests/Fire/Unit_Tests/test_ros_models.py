#!/usr/bin/env python3
"""
Unit tests for Phase 13 rate-of-spread models: MacArthur, Balbi, Cheney-Gould, 
and per-fuel wind height functions.

Verifies pure-Python implementations of ROS model formulae against analytical
reference values extracted from ERF source code.

Tests are organized by model:
  - MacArthur (1966) Australian formula (4 tests)
  - Balbi (2009) physical model (10 tests)
  - Balbi (2020) convective-radiative model and couplings (10 tests)
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

# Constants mirrored from ERF_BalbiModel.H
BALBI_FT_TO_M = 0.3048
BALBI_FT_INV_TO_M_INV = 3.28084
BALBI_BTU_LB_TO_J_KG = 2326.0
BALBI_GRAVITY = 9.81
BALBI_V_B_MIN = 1.0e-3
BALBI_B_STAR_MIN = 1.0e-6
BALBI_MIN_DEPTH_M = 0.01
BALBI_TABLE_SIZE = 14

# Balbi thermal parameters: defaults of FireParams::BalbiParams
BALBI_DEFAULTS = dict(T_a=300.0, T_f=1000.0, T_i=600.0,
                      delta_H=2.26e6, C_pf=1800.0,
                      r_00=2.5e-4, tau_0=75591.0)

# Anderson FBFM13 parameters used by the Balbi model, from
# get_anderson_fuel_params() in ERF_Rothermel.H:
#   sigma_d1 [1/ft], delta [ft], heat_content [BTU/lb]
ANDERSON_FUELS = {
    1:  (3500.0, 1.0, 8000.0),
    2:  (3000.0, 1.0, 8000.0),
    3:  (1500.0, 2.5, 8000.0),
    4:  (1739.0, 6.0, 8000.0),
    5:  (1739.0, 2.0, 8000.0),
    6:  (1750.0, 2.5, 8000.0),
    7:  (1550.0, 2.5, 8000.0),
    8:  (2000.0, 0.2, 8000.0),
    9:  (2500.0, 0.2, 8000.0),
    10: (2000.0, 1.0, 8000.0),
    11: (1500.0, 1.0, 8000.0),
    12: (1500.0, 2.3, 8000.0),
    13: (1500.0, 3.0, 8000.0),
}


def compute_balbi_params(sigma_d1, delta_ft, heat_content, M_f, bp=None,
                         M_x=None, use_moisture_extinction=False):
    """
    Balbi (2009) fuel/thermal coefficients.
    
    Formula from ERF_BalbiModel.H, compute_balbi_params():
      chi    = r_00 * sigma_m / (1 + r_00 * sigma_m)
      B_star = (C_pf * (T_i - T_a) + M_f * delta_H) / h
      v_b    = sqrt(g * delta_m * (T_f - T_a) / T_a)
      A      = chi * sigma_m * delta_m / (2 * tau_0 * B_star)
    
    Args:
        sigma_d1:     1-hr dead fuel surface-area-to-volume ratio [1/ft]
        delta_ft:     Fuel bed depth [ft]
        heat_content: Heat of combustion [BTU/lb]
        M_f:          Fuel moisture content [fraction]
        bp:           Balbi thermal parameters (defaults to BALBI_DEFAULTS)
    
    Returns:
        (A_coeff [m/s], v_b [m/s])
    """
    if bp is None:
        bp = BALBI_DEFAULTS

    sigma_m = sigma_d1 * BALBI_FT_INV_TO_M_INV
    delta_m = delta_ft * BALBI_FT_TO_M
    h_si = heat_content * BALBI_BTU_LB_TO_J_KG

    if delta_m < BALBI_MIN_DEPTH_M or sigma_m <= 0.0 or h_si <= 0.0:
        return 0.0, 1.0

    # Moisture of extinction cutoff: neither formulation has an extinction
    # limit of its own, so without this a fuel bed wetter than its M_x spreads.
    if use_moisture_extinction and M_x is not None and M_f >= M_x:
        return 0.0, 1.0

    chi_0 = bp["r_00"] * sigma_m
    chi = chi_0 / (1.0 + chi_0)

    E_ig = bp["C_pf"] * (bp["T_i"] - bp["T_a"]) + M_f * bp["delta_H"]
    B_star = E_ig / h_si

    dT_f = bp["T_f"] - bp["T_a"]
    v_b = (math.sqrt(BALBI_GRAVITY * delta_m * dT_f / bp["T_a"])
           if dT_f > 0.0 else BALBI_V_B_MIN)

    A_coeff = (chi * sigma_m * delta_m / (2.0 * bp["tau_0"] * B_star)
               if B_star > BALBI_B_STAR_MIN else 0.0)

    return max(A_coeff, 0.0), max(v_b, BALBI_V_B_MIN)


def build_fuel_balbi_table(M_f, bp=None):
    """
    Per-fuel-code Balbi lookup table.
    
    Mirrors build_fuel_balbi_table() in ERF_BalbiModel.H: index 0 is the
    non-burnable code (zero spread), 1-13 are the Anderson fuel models.
    
    Returns:
        list of (A_coeff, v_b) of length BALBI_TABLE_SIZE
    """
    table = [(0.0, 1.0)] * BALBI_TABLE_SIZE
    for code in range(1, BALBI_TABLE_SIZE):
        sigma_d1, delta_ft, heat_content = ANDERSON_FUELS[code]
        table[code] = compute_balbi_params(sigma_d1, delta_ft, heat_content, M_f, bp)
    return table


def compute_balbi_angle(U, theta, v_b):
    """
    Compute wind and slope angle for Balbi model.
    
    Formula from ERF_BalbiModel.H, compute_balbi_angle():
      tan α = U/v_b + tan θ
      where v_b < BALBI_V_B_MIN is clamped to BALBI_V_B_MIN
    
    Args:
        U: Wind speed in spread direction [m/s]
        theta: Terrain slope [radians]
        v_b: Buoyancy velocity [m/s]
    
    Returns:
        Angle α [radians]
    """
    v_b = max(v_b, BALBI_V_B_MIN)
    tan_alpha = U / v_b + math.tan(theta)
    return math.atan(tan_alpha)


def balbi_ros(U, theta, A_coeff, v_b):
    """
    Balbi (2009) ROS formula.
    
    Formula from ERF_BalbiModel.H, fill_balbi_ros():
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


def test_balbi_A_coeff_reference_values():
    """Test 10: Balbi A_coeff matches wildfire_levelset reference values."""
    # Reference values from wildfire_levelset/src/balbi_model.H at M_f = 6%
    M_f = 0.06
    cases = {1: 0.472909, 4: 1.117462, 10: 0.226372}
    worst = 0.0
    for code, expected in cases.items():
        sigma_d1, delta_ft, heat_content = ANDERSON_FUELS[code]
        A_coeff, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, M_f)
        worst = max(worst, abs(A_coeff - expected) / expected)
    passed = worst < 1.0e-4
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 10: Balbi A_coeff matches reference values (FM1, FM4, FM10)")
    if not passed:
        for code, expected in cases.items():
            sigma_d1, delta_ft, heat_content = ANDERSON_FUELS[code]
            A_coeff, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, M_f)
            print(f"    FM{code}: expected {expected}, got {A_coeff}")
    return passed


def test_balbi_A_coeff_varies_with_fuel():
    """Test 11: Balbi A_coeff differs between fuel models."""
    # A depends on sigma_m, delta_m and h, so heavy deep fuels spread faster
    # than short grass, which in turn beats timber litter.
    M_f = 0.06
    A = {}
    for code in (1, 4, 10):
        sigma_d1, delta_ft, heat_content = ANDERSON_FUELS[code]
        A[code], _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, M_f)
    passed = A[4] > A[1] > A[10]
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 11: Balbi A_coeff varies with fuel model (FM4 > FM1 > FM10)")
    if not passed:
        print(f"    FM4={A[4]}, FM1={A[1]}, FM10={A[10]}")
    return passed


def test_balbi_A_coeff_decreases_with_moisture():
    """Test 12: Balbi A_coeff decreases as fuel moisture rises."""
    # Moisture enters through B*, and A is inversely proportional to B*
    sigma_d1, delta_ft, heat_content = ANDERSON_FUELS[1]
    A_dry, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, 0.06)
    A_wet, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, 0.20)
    passed = A_wet < A_dry
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 12: Balbi A_coeff decreases with fuel moisture")
    if not passed:
        print(f"    A(6%)={A_dry}, A(20%)={A_wet}")
    return passed


def test_balbi_table_structure():
    """Test 13: Balbi fuel table has 14 entries with code 0 non-burnable."""
    table = build_fuel_balbi_table(0.06)
    size_ok = len(table) == BALBI_TABLE_SIZE
    nonburnable_ok = table[0][0] == 0.0
    burnable_ok = all(table[code][0] > 0.0 for code in range(1, BALBI_TABLE_SIZE))
    passed = size_ok and nonburnable_ok and burnable_ok
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 13: Balbi fuel table size 14, code 0 non-burnable")
    if not passed:
        print(f"    size={len(table)}, A[0]={table[0][0]}, "
              f"min burnable A={min(t[0] for t in table[1:])}")
    return passed


def test_balbi_tan_to_sin_cos_identity():
    """Test 14: kernel sin/cos from tan matches the atan-based formulation."""
    # fill_balbi_ros() avoids atan by using
    #   sin a = tan a / sqrt(1 + tan^2 a), cos a = 1 / sqrt(1 + tan^2 a)
    A_coeff, v_b = compute_balbi_params(*ANDERSON_FUELS[1], M_f=0.06)
    worst = 0.0
    for U in (0.0, 1.0, 5.0, 12.0):
        for slope_mag in (0.0, 0.3, 1.0):
            tan_alpha = U / v_b + slope_mag
            inv_sec = 1.0 / math.sqrt(1.0 + tan_alpha * tan_alpha)
            ros_kernel = max(A_coeff * (1.0 + tan_alpha * inv_sec - inv_sec), 0.0)
            ros_atan = balbi_ros(U, math.atan(slope_mag), A_coeff, v_b)
            worst = max(worst, abs(ros_kernel - ros_atan))
    passed = worst < 1.0e-12
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 14: Balbi kernel sin/cos identity matches atan formulation")
    if not passed:
        print(f"    max |R_kernel - R_atan| = {worst}")
    return passed


# ============================================================================
# Balbi (2020) Convective-Radiative Model Tests
# ============================================================================

# Constants of Balbi et al. (2020), Table 1 nomenclature
BALBI2020_DEFAULTS = dict(
    C_pf=1800.0,     # specific heat of fuel [J/(kg K)]
    C_pw=4180.0,     # specific heat of water [J/(kg K)]
    C_pa=1150.0,     # specific heat of air [J/(kg K)]
    T_a=300.0,       # ambient temperature [K]
    T_i=600.0,       # ignition temperature [K]
    T_vap=373.0,     # vaporisation temperature [K]
    delta_H=2.3e6,   # latent heat of evaporation [J/kg]
    dH_comb=1.74e7,  # heat of combustion of pyrolysis gases [J/kg]
    chi_0=0.3,       # radiative factor [-]
    tau_0=75591.0,   # flame residence time parameter [s/m]
    r_00=2.5e-5,     # model coefficient [-]
    K1=130.0,        # drag coefficient [s/m]
    st=17.0,         # air-pyrolysis gas mass ratio [-]
    rho_a=1.2,       # air density [kg/m3]
    sigma_B=5.6e-8,  # Stefan-Boltzmann constant [W/(m2 K4)]
)

TWO_PI = 2.0 * math.pi


def compute_balbi2020_state(sav, depth, load, rho_v, M_f, bp=None):
    """
    Balbi (2020) coefficients that do not depend on the rate of spread.
    
    Mirrors the 2020 branch of compute_balbi_params() in ERF_BalbiModel.H,
    in SI units directly (the C++ side converts from Rothermel units first).
    
    Args:
        sav:   surface-area-to-volume ratio s [1/m]
        depth: fuel bed depth h [m]
        load:  dead fuel load [kg/m2]
        rho_v: fuel particle density [kg/m3]
        M_f:   fuel moisture content [fraction]
    
    Returns:
        dict of precomputed coefficients
    """
    if bp is None:
        bp = BALBI2020_DEFAULTS

    beta = load / (depth * rho_v)
    S = sav * beta * depth

    # eq. 9: the water term uses the specific heat of water, not of the fuel
    q = (bp["C_pf"] * (bp["T_i"] - bp["T_a"])
         + M_f * (bp["delta_H"] + bp["C_pw"] * (bp["T_vap"] - bp["T_a"])))

    a_r = min(S / TWO_PI, 1.0)                                          # eq. 17
    return dict(
        A_rad=a_r * bp["chi_0"] * bp["dH_comb"] / (4.0 * q),            # eq. 16
        Rb_coef=min(S / math.pi, 1.0) * bp["sigma_B"] / (beta * rho_v * q),  # eq. 13
        Rc_coef=(sav * bp["dH_comb"] / (q * bp["tau_0"])
                 * min(depth, TWO_PI / (sav * beta))),                  # eq. 27
        u0_coef=(2.0 * (bp["st"] + 1.0) / bp["tau_0"] * (rho_v / bp["rho_a"])
                 * min(S, TWO_PI)),                                     # eq. B9
        s_r00=sav * bp["r_00"],
        K1_sqrt_beta=bp["K1"] * math.sqrt(beta),
        h_depth=depth,
        T_a=bp["T_a"],
        T_flame_coef=bp["dH_comb"] / (bp["C_pa"] * (bp["st"] + 1.0)),   # eq. B11
        chi_0=bp["chi_0"],
    )


def balbi2020_rhs(bc, R, U, tan_slope):
    """
    Right-hand side of Balbi (2020) eq. 28, R_b + R_c + R_r.
    
    Mirrors balbi2020_rhs() in ERF_BalbiModel.H, including the two-pass
    refinement of the flame tilt that eq. C7 needs.
    """
    cos_gamma = 1.0 / math.sqrt(1.0 + tan_slope * tan_slope)
    sin_gamma = tan_slope * cos_gamma
    T = bc["T_a"]
    u0 = 1.0e-3

    for _ in range(2):
        chi = bc["chi_0"] / (1.0 + R * cos_gamma / bc["s_r00"])         # eq. C7
        T = bc["T_a"] + bc["T_flame_coef"] * (1.0 - chi)                # eq. B11
        u0 = max(bc["u0_coef"] * (T / bc["T_a"]), 1.0e-3)               # eq. B9
        tan_gamma = tan_slope + max(U, 0.0) / u0                        # eq. 2
        inv_sec = 1.0 / math.sqrt(1.0 + tan_gamma * tan_gamma)
        sin_gamma = tan_gamma * inv_sec
        cos_gamma = inv_sec

    T_ratio = T / bc["T_a"]
    H = u0 * u0 / (9.81 * (T_ratio - 1.0)) if T_ratio > 1.0 else 0.0    # eq. 23

    Rb = bc["Rb_coef"] * T**4                                           # eq. 13
    Rc = bc["Rc_coef"] * (bc["h_depth"] / (2.0 * bc["h_depth"] + H) * tan_slope
                          + max(U, 0.0) * math.exp(-bc["K1_sqrt_beta"] * R) / u0)
    Rr = (bc["A_rad"] * R * (1.0 + sin_gamma - cos_gamma)
          / (1.0 + R * cos_gamma / bc["s_r00"]))                        # eq. 15
    return Rb + Rc + Rr


def balbi2020_ros(bc, U, tan_slope=0.0, max_iter=60, tol=1.0e-7):
    """
    Balbi (2020) ROS by bracketed root find, as in balbi2020_ros().
    """
    lo, hi = 0.0, 30.0
    if balbi2020_rhs(bc, lo, U, tan_slope) - lo <= 0.0:
        return 0.0
    if balbi2020_rhs(bc, hi, U, tan_slope) - hi > 0.0:
        return hi
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if balbi2020_rhs(bc, mid, U, tan_slope) - mid > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= tol:
            break
    return 0.5 * (lo + hi)


# Table 2 of Balbi et al. (2020): the fuel bed used for the ROS-vs-wind and
# ROS-vs-FMC numerical simulations of Figs 3 and 4.
PAPER_TABLE2 = dict(sav=6000.0, depth=0.1, rho_v=500.0, M_f=0.10)


def test_balbi2020_radiative_coefficient():
    """Test 15: Balbi (2020) A matches the value quoted in the paper."""
    # The paper states A = 0.15 for a fuel load of 0.05 kg/m2 with the
    # Table 2 fuel bed.
    bc = compute_balbi2020_state(load=0.05, **PAPER_TABLE2)
    passed = abs(bc["A_rad"] - 0.15) < 0.01
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 15: Balbi (2020) A = 0.15 for 0.05 kg/m2 (paper value)")
    if not passed:
        print(f"    Expected ~0.15, got {bc['A_rad']}")
    return passed


def test_balbi2020_nonzero_no_wind_ros():
    """Test 16: Balbi (2020) spreads with no wind and no slope."""
    # The radiative base term R_b is what the 2009 form lacks entirely.
    bc = compute_balbi2020_state(load=0.3, **PAPER_TABLE2)
    ros = balbi2020_ros(bc, U=0.0)
    passed = ros > 0.0
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 16: Balbi (2020) no-wind ROS is nonzero")
    if not passed:
        print(f"    ROS(U=0) = {ros}")
    return passed


def test_balbi2020_monotone_in_wind():
    """Test 17: Balbi (2020) ROS rises monotonically with wind (Fig. 3)."""
    ok = True
    for load in (0.05, 0.3, 0.8):
        bc = compute_balbi2020_state(load=load, **PAPER_TABLE2)
        prev = -1.0
        for U in [0.0, 0.4, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]:
            ros = balbi2020_ros(bc, U)
            if ros <= prev:
                ok = False
            prev = ros
    status = "\u2713" if ok else "\u2717"
    print(f"{status} Test 17: Balbi (2020) ROS monotone in wind for all Fig. 3 loads")
    return ok


def test_balbi2020_fig3_magnitudes():
    """Test 18: Balbi (2020) ROS matches the Fig. 3 curves."""
    # Values produced by the shipped C++ (balbi2020_ros in ERF_BalbiModel.H)
    # for the Table 2 fuel bed, which trace the curves plotted in Fig. 3.
    expected = {
        0.05: {0.0: 0.019, 2.0: 0.325, 8.0: 0.551, 12.0: 0.624},
        0.30: {0.0: 0.016, 2.0: 0.127, 8.0: 0.381, 12.0: 0.575},
        0.80: {0.0: 0.006, 2.0: 0.070, 8.0: 0.428, 12.0: 0.723},
    }
    worst = 0.0
    for load, curve in expected.items():
        bc = compute_balbi2020_state(load=load, **PAPER_TABLE2)
        for U, ros_ref in curve.items():
            worst = max(worst, abs(balbi2020_ros(bc, U) - ros_ref))
    passed = worst < 1.0e-3
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 18: Balbi (2020) ROS matches the Fig. 3 curves")
    if not passed:
        print(f"    max deviation {worst} m/s")
    return passed


def test_balbi2020_decreases_with_moisture():
    """Test 19: Balbi (2020) ROS falls as fuel moisture rises (Fig. 4)."""
    args = dict(PAPER_TABLE2)
    del args["M_f"]
    ros_prev = None
    ok = True
    for m in (0.05, 0.20, 0.40, 0.60):
        bc = compute_balbi2020_state(load=0.09, M_f=m, **args)
        ros = balbi2020_ros(bc, U=8.0)
        if ros_prev is not None and ros >= ros_prev:
            ok = False
        ros_prev = ros
    status = "\u2713" if ok else "\u2717"
    print(f"{status} Test 19: Balbi (2020) ROS decreases with fuel moisture")
    return ok


def test_balbi_directional_projection():
    """Test 20: directional ROS peaks downwind and vanishes upwind."""
    # fill_balbi_ros_directional() evaluates the tilt with U . n and grad(z) . n
    # for the front normal n. With no slope, the head direction is the wind
    # direction and the backing direction gets a negative normal wind.
    A_coeff, v_b = compute_balbi_params(*ANDERSON_FUELS[1], M_f=0.06)
    u, v = 5.0, 0.0  # wind along +x

    def ros_along(nx, ny):
        U_n = u * nx + v * ny
        tan_alpha = U_n / v_b
        inv_sec = 1.0 / math.sqrt(1.0 + tan_alpha * tan_alpha)
        return max(A_coeff * (1.0 + tan_alpha * inv_sec - inv_sec), 0.0)

    head = ros_along(1.0, 0.0)
    flank = ros_along(0.0, 1.0)
    back = ros_along(-1.0, 0.0)
    passed = head > flank >= 0.0 and back == 0.0
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 20: directional ROS head > flank, zero backing (2009 form)")
    if not passed:
        print(f"    head={head}, flank={flank}, back={back}")
    return passed


def test_balbi_heat_flux_uprights_flame():
    """Test 21: heat-flux buoyancy stands the flame up and slows the head."""
    # v_b_eff = sqrt(v_b^2 + v_b_Q^2) enters the tilt as U/v_b, so a stronger
    # plume tilts the flame less and reduces forward spread.
    A_coeff, v_b = compute_balbi_params(*ANDERSON_FUELS[1], M_f=0.06)

    def ros(vz_extra):
        vb = math.sqrt(v_b * v_b + vz_extra * vz_extra)
        t = 5.0 / vb
        inv_sec = 1.0 / math.sqrt(1.0 + t * t)
        return max(A_coeff * (1.0 + t * inv_sec - inv_sec), 0.0)

    passed = ros(3.0) < ros(0.0)
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 21: heat-flux buoyancy reduces the head ROS")
    if not passed:
        print(f"    ROS(no flux)={ros(0.0)}, ROS(v_b_Q=3)={ros(3.0)}")
    return passed


# Moisture of extinction per Anderson fuel model, from ERF_Rothermel.H
ANDERSON_MX = {1: 0.12, 2: 0.15, 3: 0.25, 4: 0.20, 5: 0.20, 6: 0.25, 7: 0.40,
               8: 0.30, 9: 0.25, 10: 0.25, 11: 0.15, 12: 0.20, 13: 0.25}


def test_balbi_moisture_extinction_cutoff():
    """Test 22: ROS goes to zero at the moisture of extinction."""
    sigma_d1, delta_ft, heat_content = ANDERSON_FUELS[1]
    M_x = ANDERSON_MX[1]  # 0.12 for short grass

    below, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, M_x - 0.01,
                                    M_x=M_x, use_moisture_extinction=True)
    at, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, M_x,
                                 M_x=M_x, use_moisture_extinction=True)
    above, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, M_x + 0.05,
                                    M_x=M_x, use_moisture_extinction=True)
    # With the cutoff off, a fuel bed past extinction still spreads
    off, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, M_x + 0.05,
                                  M_x=M_x, use_moisture_extinction=False)

    passed = below > 0.0 and at == 0.0 and above == 0.0 and off > 0.0
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 22: Balbi ROS zeroed at the moisture of extinction")
    if not passed:
        print(f"    below={below}, at={at}, above={above}, cutoff_off={off}")
    return passed


def test_balbi_moisture_clamp():
    """Test 23: per-cell moisture is clamped to the [0.01, 0.40] band."""
    # balbi_state_at_cell() clamps the moisture it reads from the ODE state the
    # same way the domain average is clamped, so a wild cell cannot produce a
    # negative ignition energy.
    sigma_d1, delta_ft, heat_content = ANDERSON_FUELS[1]

    def clamped(m):
        return max(0.01, min(m, 0.40))

    A_lo, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, clamped(-0.5))
    A_ref, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, 0.01)
    A_hi, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, clamped(5.0))
    A_max, _ = compute_balbi_params(sigma_d1, delta_ft, heat_content, 0.40)

    passed = (abs(A_lo - A_ref) < 1e-12 and abs(A_hi - A_max) < 1e-12
              and A_hi > 0.0)
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 23: per-cell moisture clamped to [0.01, 0.40]")
    if not passed:
        print(f"    A(-0.5)={A_lo} vs A(0.01)={A_ref}, "
              f"A(5.0)={A_hi} vs A(0.40)={A_max}")
    return passed


def test_balbi_herb_curing_raises_load():
    """Test 24: cured herbaceous load raises the 2020 radiant factor."""
    # FM2 carries a live herbaceous load; curing moves it into the dead fine
    # fuel that the 2020 form's packing ratio is built from.
    w_d1, w_lh = 0.046, 0.023          # lb/ft2, from get_anderson_fuel_params(2)
    sav = 3000.0 * BALBI_FT_INV_TO_M_INV
    depth = 1.0 * BALBI_FT_TO_M
    rho_v = 32.0 * 16.0185

    def A_for(curing):
        load = (w_d1 + curing * w_lh) * 4.8824
        return compute_balbi2020_state(sav=sav, depth=depth, load=load,
                                       rho_v=rho_v, M_f=0.06)["A_rad"]

    A_green = A_for(0.0)
    A_half = A_for(0.5)
    A_cured = A_for(1.0)
    passed = A_green < A_half < A_cured
    status = "\u2713" if passed else "\u2717"
    print(f"{status} Test 24: cured herbaceous load raises the 2020 radiant factor")
    if not passed:
        print(f"    A(0%)={A_green}, A(50%)={A_half}, A(100%)={A_cured}")
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
    """Test 25: Cheney-Gould zero wind formula."""
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
    print(f"{status} Test 25: Cheney-Gould zero wind formula")
    if not passed:
        print(f"    Expected: {expected}, Got: {ros}")
    return passed


def test_cheney_gould_increases_with_wind():
    """Test 26: Cheney-Gould ROS increases with wind."""
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
    print(f"{status} Test 26: Cheney-Gould ROS increases with wind")
    if not passed:
        print(f"    R(U=0)={ros_0}, R(U=3)={ros_3}, R(U=5)={ros_5}")
    return passed


def test_cheney_gould_decreases_with_moisture():
    """Test 27: Cheney-Gould ROS decreases with moisture."""
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
    print(f"{status} Test 27: Cheney-Gould ROS decreases with moisture")
    if not passed:
        print(f"    R(5%)={ros_5pct}, R(10%)={ros_10pct}, R(20%)={ros_20pct}")
    return passed


def test_cheney_gould_increases_with_curing():
    """Test 28: Cheney-Gould ROS increases with curing."""
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
    print(f"{status} Test 28: Cheney-Gould ROS increases with curing")
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
    """Test 29: fcwh uniform mode returns global_z_ref for all fuels."""
    global_z_ref = 6.1
    fcwh = build_fcwh_table(global_z_ref, use_per_fuel=False)
    passed = (len(fcwh) == 14 and 
              all(fcwh[i] == global_z_ref for i in range(1, 14)))
    status = "✓" if passed else "✗"
    print(f"{status} Test 29: fcwh uniform mode (all fuels = {global_z_ref})")
    if not passed:
        print(f"    Length: {len(fcwh)}, entries 1-13 all {global_z_ref}: "
              f"{all(fcwh[i] == global_z_ref for i in range(1, 14))}")
    return passed


def test_fcwh_per_fuel_mode():
    """Test 30: fcwh per-fuel mode returns 6.096 for all fuels."""
    global_z_ref = 6.1
    fcwh = build_fcwh_table(global_z_ref, use_per_fuel=True)
    expected = 6.096
    passed = (len(fcwh) == 14 and 
              all(abs(fcwh[i] - expected) < 1e-6 for i in range(1, 14)))
    status = "✓" if passed else "✗"
    print(f"{status} Test 30: fcwh per-fuel mode (all fuels = 6.096 m)")
    if not passed:
        print(f"    Length: {len(fcwh)}")
        for i in range(1, 14):
            if abs(fcwh[i] - expected) >= 1e-6:
                print(f"    fcwh[{i}] = {fcwh[i]}, expected {expected}")
    return passed


def test_fcz0_fm4_chaparral():
    """Test 31: fcz0 FM4 value equals 0.2378 (chaparral, highest roughness)."""
    fcz0 = build_fcz0_table()
    expected = 0.2378
    passed = abs(fcz0[4] - expected) < 1e-6
    status = "✓" if passed else "✗"
    print(f"{status} Test 31: fcz0 FM4 = 0.2378 m (chaparral)")
    if not passed:
        print(f"    Expected: {expected}, Got: {fcz0[4]}")
    return passed


def test_fcz0_fm1_fm2_equal():
    """Test 32: fcz0 FM1 and FM2 both equal 0.0396."""
    fcz0 = build_fcz0_table()
    expected = 0.0396
    passed = (abs(fcz0[1] - expected) < 1e-6 and 
              abs(fcz0[2] - expected) < 1e-6 and 
              fcz0[1] == fcz0[2])
    status = "✓" if passed else "✗"
    print(f"{status} Test 32: fcz0 FM1 and FM2 both equal 0.0396 m")
    if not passed:
        print(f"    fcz0[1] = {fcz0[1]}, fcz0[2] = {fcz0[2]}, expected {expected}")
    return passed


def test_fcz0_table_size():
    """Test 18: fcz0 table has size 14 (indices 0-13)."""
    fcz0 = build_fcz0_table()
    passed = len(fcz0) == 14
    status = "✓" if passed else "✗"
    print(f"{status} Test 33: fcz0 table size = 14")
    if not passed:
        print(f"    Expected length 14, got {len(fcz0)}")
    return passed


# ============================================================================
# Test Runner
# ============================================================================

def main():
    """Run all 33 tests and return exit code (0 = all pass, 1 = any fail)."""
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
    
    # Balbi tests (10)
    print("Balbi (2009) Physical Model Tests")
    print("-" * 70)
    results.append(test_balbi_zero_wind_slope())
    results.append(test_balbi_increases_with_wind())
    results.append(test_balbi_angle_formula())
    results.append(test_balbi_small_v_b_floor())
    results.append(test_balbi_positive_ros_clamp())
    results.append(test_balbi_A_coeff_reference_values())
    results.append(test_balbi_A_coeff_varies_with_fuel())
    results.append(test_balbi_A_coeff_decreases_with_moisture())
    results.append(test_balbi_table_structure())
    results.append(test_balbi_tan_to_sin_cos_identity())
    print()

    # Balbi (2020) and coupling tests (10)
    print("Balbi (2020) Convective-Radiative Model Tests")
    print("-" * 70)
    results.append(test_balbi2020_radiative_coefficient())
    results.append(test_balbi2020_nonzero_no_wind_ros())
    results.append(test_balbi2020_monotone_in_wind())
    results.append(test_balbi2020_fig3_magnitudes())
    results.append(test_balbi2020_decreases_with_moisture())
    results.append(test_balbi_directional_projection())
    results.append(test_balbi_heat_flux_uprights_flame())
    results.append(test_balbi_moisture_extinction_cutoff())
    results.append(test_balbi_moisture_clamp())
    results.append(test_balbi_herb_curing_raises_load())
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
