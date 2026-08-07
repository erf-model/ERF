#!/usr/bin/env python3
"""
Phase 6 Two-Stream Radiation: Time-Integration Timing Check

Checks:
  1) Expected diagnostic row count for configured cadence
  2) Per-step multiplicity (rows per step)
  3) SW_TOA matches S0*cos(zenith)
  4) heating_rate_max finite and nonzero
  5) heating_rate_max stability (coefficient of variation)

This script is designed for current Phase 6 behavior where diagnostics may be
written multiple times per coarse step (e.g., pre/post dycore).
"""

import csv
import math
import statistics
import sys
from collections import Counter
from pathlib import Path

# ---------------------------------------------------------------------
# User-tunable test configuration
# ---------------------------------------------------------------------
CSV_PATH = Path("radiation_phase6_timing_diag.dat")

S0 = 1361.0
SOLAR_ZENITH_DEG = 60.0

DT = 0.5
STOP_TIME = 5.0

# Phase 6 current expected behavior: pre + post dycore diagnostic append
DIAG_CALLS_PER_STEP = 2

# Tolerances
ROW_TOLERANCE = DIAG_CALLS_PER_STEP      # allow small startup/teardown variation
SW_TOA_ABS_TOL = 1.0e-6
CV_STABILITY_TOL = 0.05                  # 5% variation threshold
HEATING_NONZERO_TOL = 1.0e-12


def fail(msg: str):
    print(f"[FAIL] {msg}")
    raise SystemExit(1)


def pass_line(msg: str):
    print(f"  {msg} [PASS]")


def read_rows(csv_path: Path):
    if not csv_path.exists():
        fail(f"Diagnostic CSV not found: {csv_path}")

    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        fail("Diagnostic CSV is empty.")

    required_cols = [
        "step",
        "time",
        "SW_TOA",
        "heating_rate_max",
    ]
    missing = [c for c in required_cols if c not in rows[0]]
    if missing:
        fail(f"Missing required columns in CSV: {missing}")

    return rows


def to_float(row, key):
    try:
        return float(row[key])
    except Exception:
        fail(f"Could not parse float column '{key}' from row: {row}")


def to_int(row, key):
    try:
        return int(float(row[key]))
    except Exception:
        fail(f"Could not parse int column '{key}' from row: {row}")


def main():
    print("\n======================================================================")
    print("Phase 6 Two-Stream Radiation: Time-Integration Timing Check")
    print("======================================================================\n")

    cosz = math.cos(math.radians(SOLAR_ZENITH_DEG))
    expected_toa = S0 * cosz
    nsteps = int(round(STOP_TIME / DT))
    expected_rows = DIAG_CALLS_PER_STEP * nsteps

    print("Test Parameters:")
    print(f"  Solar constant S0 = {S0:.2f} W/m^2")
    print(f"  Solar zenith angle = {SOLAR_ZENITH_DEG:.1f}°")
    print(f"  cos(zenith) = {cosz:.4f}")
    print(f"  Expected TOA flux = {expected_toa:.4f} W/m^2")
    print(f"  Fixed dt = {DT} s, stop_time = {STOP_TIME} s (expect ~{nsteps} coarse steps)")
    print(f"  Diagnostic calls per step = {DIAG_CALLS_PER_STEP} (expect ~{expected_rows} rows)\n")

    rows = read_rows(CSV_PATH)
    nrows = len(rows)
    print(f"Diagnostic CSV rows found: {nrows}")

    # -----------------------------------------------------------------
    # 1) Row count check
    # -----------------------------------------------------------------
    if abs(nrows - expected_rows) <= ROW_TOLERANCE:
        pass_line(f"Row count check: {nrows} rows (expected ~{expected_rows})")
    else:
        fail(
            f"Expected approximately {expected_rows} diagnostic rows "
            f"(±{ROW_TOLERANCE}), but found {nrows}."
        )

    # -----------------------------------------------------------------
    # 2) Per-step multiplicity check
    # -----------------------------------------------------------------
    steps = [to_int(r, "step") for r in rows]
    step_counts = Counter(steps)

    bad_counts = {s: c for s, c in sorted(step_counts.items()) if c != DIAG_CALLS_PER_STEP}
    if not bad_counts:
        pass_line(f"Step multiplicity check ({DIAG_CALLS_PER_STEP} rows/step)")
    else:
        fail(
            "Unexpected diagnostic multiplicity by step: "
            + ", ".join(f"step {s}: {c}" for s, c in bad_counts.items())
        )

    # -----------------------------------------------------------------
    # 3) SW_TOA accuracy check
    # -----------------------------------------------------------------
    sw_toa_vals = [to_float(r, "SW_TOA") for r in rows]
    sw_bad = [
        (i, v) for i, v in enumerate(sw_toa_vals)
        if abs(v - expected_toa) > SW_TOA_ABS_TOL
    ]
    if not sw_bad:
        pass_line(f"SW_TOA accuracy check (all {nrows} rows)")
    else:
        i0, v0 = sw_bad[0]
        fail(
            f"SW_TOA mismatch at row {i0}: got {v0:.12g}, "
            f"expected {expected_toa:.12g}, tol={SW_TOA_ABS_TOL}"
        )

    # -----------------------------------------------------------------
    # 4) heating_rate_max finite and nonzero
    # -----------------------------------------------------------------
    hvals = [to_float(r, "heating_rate_max") for r in rows]

    if all(math.isfinite(v) for v in hvals):
        pass_line("heating_rate_max finite (no NaN/Inf) check")
    else:
        fail("Found non-finite heating_rate_max (NaN or Inf).")

    mean_h = statistics.mean(hvals)
    final_h = hvals[-1]

    if any(abs(v) > HEATING_NONZERO_TOL for v in hvals):
        pass_line(
            "heating_rate_max nonzero check\n"
            f"    (mean heating_rate_max = {mean_h:.6e} K/s, final = {final_h:.6e} K/s)"
        )
    else:
        fail("All heating_rate_max values are effectively zero.")

    # -----------------------------------------------------------------
    # 5) Stability check via coefficient of variation
    # -----------------------------------------------------------------
    if len(hvals) > 1:
        stdev_h = statistics.pstdev(hvals)
        cv = abs(stdev_h / mean_h) if abs(mean_h) > 0 else float("inf")
    else:
        cv = 0.0

    if cv <= CV_STABILITY_TOL:
        pass_line(f"heating_rate_max stability (CV = {cv:.4f})")
    else:
        fail(
            f"heating_rate_max is too variable: CV={cv:.6f}, "
            f"threshold={CV_STABILITY_TOL:.6f}"
        )

    # -----------------------------------------------------------------
    # Optional informational check: duplicate (step,time) pairs
    # (May be expected if cadence includes pre/post at same step.)
    # -----------------------------------------------------------------
    step_time = [(to_int(r, "step"), to_float(r, "time")) for r in rows]
    uniq_step_time = len(set(step_time))
    if uniq_step_time == len(step_time):
        pass_line("Unique (step,time) pairs check")
    else:
        print(
            "  Unique (step,time) pairs check [INFO]\n"
            "    Duplicate (step,time) pairs detected. "
            "This can be valid when multiple diagnostic call-sites are used per step."
        )

    print("\n======================================================================")
    print("TEST PASSED")
    print("======================================================================\n")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        print("\n======================================================================")
        print("TEST FAILED")
        print("======================================================================\n")
        raise
