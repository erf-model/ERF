#!/usr/bin/env python3
"""
Phase 8 Benchmark Suite: Central Tolerance Configuration

This module defines all tolerances used across the benchmark suite to avoid
magic numbers scattered throughout the codebase. All metric checks should
reference tolerances defined here.
"""

# =========================================================================
# FLUX TOLERANCES
# =========================================================================

# Shortwave TOA flux relative error tolerance (%)
SW_TOA_RELATIVE_TOL_PCT = 0.1

# Shortwave surface flux relative error tolerance (%)
SW_SURFACE_RELATIVE_TOL_PCT = 1.0

# Longwave surface net flux relative error tolerance (%)
LW_NET_SURFACE_RELATIVE_TOL_PCT = 1.0

# Absolute tolerance for SW/LW flux comparisons (W/m^2)
# Used when relative error doesn't make sense (e.g., very small fluxes)
FLUX_ABSOLUTE_TOL_W_M2 = 1.0e-6

# =========================================================================
# HEATING RATE TOLERANCES
# =========================================================================

# Coefficient of variation threshold for heating_rate_max stability (%)
# Checks that max heating rate doesn't oscillate too much across steps
HEATING_CV_UPPER_BOUND = 0.05  # 5%

# Absolute threshold to consider heating rate nonzero (K/s)
HEATING_NONZERO_TOL = 1.0e-12

# =========================================================================
# CADENCE TOLERANCES
# =========================================================================

# Row count tolerance (absolute number of rows)
# Allows small variation due to startup/teardown effects
ROW_COUNT_ABS_TOL = 2

# Acceptable range for rows per step
# "both" mode: expect 2 (pre + post)
# "pre_only" or "post_only": expect 1
# Actual value depends on diag_callsite_mode configured for each case
ROWS_PER_STEP_TOL = 1  # Allow ±1 row variation

# =========================================================================
# MONOTONICITY/STABILITY TOLERANCES
# =========================================================================

# Allow small oscillations in metrics (e.g., due to numerical precision)
# Used for checks like "flux should be bounded"
STABILITY_RELATIVE_TOL = 0.01  # 1%

# =========================================================================
# HELPER FUNCTIONS
# =========================================================================

def check_relative_error(actual, expected, tol_pct):
    """
    Check if actual value is within relative tolerance of expected.
    
    Args:
        actual: Computed/observed value
        expected: Reference/analytical value
        tol_pct: Tolerance as percentage (e.g., 1.0 for 1%)
    
    Returns:
        Tuple (is_pass, error_pct)
    """
    if expected == 0:
        rel_error = abs(actual - expected)
    else:
        rel_error = abs((actual - expected) / expected) * 100
    
    return (rel_error <= tol_pct, rel_error)


def check_absolute_error(actual, expected, tol_abs):
    """
    Check if actual value is within absolute tolerance of expected.
    
    Args:
        actual: Computed/observed value
        expected: Reference/analytical value
        tol_abs: Absolute tolerance
    
    Returns:
        Tuple (is_pass, error_abs)
    """
    error = abs(actual - expected)
    return (error <= tol_abs, error)


def describe_tolerance(name, value, unit=""):
    """Helper to describe a tolerance setting."""
    return f"{name}: {value} {unit}".strip()


# Print summary of tolerances
def print_tolerance_summary():
    """Print a summary of all configured tolerances."""
    print("\n" + "=" * 70)
    print("BENCHMARK SUITE TOLERANCE CONFIGURATION")
    print("=" * 70)
    print("\nFLUX TOLERANCES:")
    print(f"  SW_TOA relative error: {SW_TOA_RELATIVE_TOL_PCT}%")
    print(f"  SW_surface relative error: {SW_SURFACE_RELATIVE_TOL_PCT}%")
    print(f"  LW_net_surface relative error: {LW_NET_SURFACE_RELATIVE_TOL_PCT}%")
    print(f"  Flux absolute: {FLUX_ABSOLUTE_TOL_W_M2} W/m^2")
    print("\nHEATING RATE TOLERANCES:")
    print(f"  heating_rate_max CV upper bound: {HEATING_CV_UPPER_BOUND * 100}%")
    print(f"  heating_rate_max nonzero threshold: {HEATING_NONZERO_TOL} K/s")
    print("\nCADENCE TOLERANCES:")
    print(f"  Row count absolute: ±{ROW_COUNT_ABS_TOL}")
    print(f"  Rows per step: ±{ROWS_PER_STEP_TOL}")
    print("\nSTABILITY TOLERANCES:")
    print(f"  Relative stability tolerance: {STABILITY_RELATIVE_TOL * 100}%")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    print_tolerance_summary()
