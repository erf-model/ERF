#!/usr/bin/env python3
"""
Phase 7 Two-Stream Radiation: Time-Integration Timing Check with Diagnostics Controls

Checks:
  1) Expected diagnostic row count for configured cadence and diagnostics mode
  2) Per-step multiplicity (rows per step)
  3) SW_TOA matches S0*cos(zenith)
  4) heating_rate_max finite and nonzero
  5) heating_rate_max stability (coefficient of variation)
  6) Call-site filtering behavior based on diag_callsite_mode

This script is enhanced for Phase 7 diagnostics controls, accounting for:
  - diag_enable: master switch (if false, may be no diagnostics file)
  - diag_callsite_mode: "both" (pre+post), "pre_only", or "post_only"
  - Other output stream controls (stdout, tagged, regtest, csv)

When diag_callsite_mode is configured, row count expectations change:
  - "both": ~2 rows per step (pre + post)
  - "pre_only": ~1 row per step (pre only)
  - "post_only": ~1 row per step (post only)
"""

import csv
import math
import statistics
import sys
from collections import Counter
from pathlib import Path

# ---------------------------------------------------------------------
# User-tunable test configuration
# Use environment variables or hardcoded defaults
# ---------------------------------------------------------------------
CSV_PATH = Path("radiation_phase6_timing_diag.dat")

S0 = 1361.0
SOLAR_ZENITH_DEG = 60.0

DT = 0.5
STOP_TIME = 5.0

# Phase 7 diagnostics configuration
# Can be set via environment or command line
DIAG_CALLSITE_MODE = "both"  # "both", "pre_only", or "post_only"
DIAG_ENABLE = True

# Compute expected diagnostic calls per step based on mode
if not DIAG_ENABLE:
    DIAG_CALLS_PER_STEP = 0
elif DIAG_CALLSITE_MODE == "both":
    DIAG_CALLS_PER_STEP = 2  # pre + post
else:  # "pre_only" or "post_only"
    DIAG_CALLS_PER_STEP = 1

# Tolerances
ROW_TOLERANCE = max(2, DIAG_CALLS_PER_STEP)  # allow small startup/teardown variation
SW_TOA_ABS_TOL = 1.0e-6
CV_STABILITY_TOL = 0.05                      # 5% variation threshold
HEATING_NONZERO_TOL = 1.0e-12


def fail(msg: str):
    print(f"[FAIL] {msg}")
    raise SystemExit(1)


def pass_line(msg: str):
    print(f"  {msg} [PASS]")


def read_rows(csv_path: Path):
    if not csv_path.exists():
        if not DIAG_ENABLE:
            print("  No diagnostic file found (expected for diag_enable=false)")
            return []
        fail(f"Diagnostic CSV not found: {csv_path}")

    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows and DIAG_ENABLE:
        fail("Diagnostic CSV is empty but diag_enable=true.")

    required_cols = [
        "step",
        "time",
        "call_site",
        "SW_TOA",
        "heating_rate_max",
    ]
    if rows:
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
    print("Phase 7 Two-Stream Radiation: Time-Integration Timing Check")
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
    print(f"  Diagnostics enabled: {DIAG_ENABLE}")
    print(f"  Call-site mode: {DIAG_CALLSITE_MODE}")
    print(f"  Diagnostic calls per step = {DIAG_CALLS_PER_STEP} (expect ~{expected_rows} rows)\n")

    rows = read_rows(CSV_PATH)

    # Handle case where diagnostics are disabled
    if not DIAG_ENABLE:
        if not rows:
            pass_line("No diagnostic rows emitted (diag_enable=false)")
        else:
            print("  [INFO] Diagnostic rows present despite diag_enable=false (may be artifact)")
        print("\n======================================================================")
        print("TEST PASSED")
        print("======================================================================\n")
        return

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
    # 3) Call-site mode validation
    # -----------------------------------------------------------------
    call_sites = [r.get("call_site", "") for r in rows]
    has_pre = any("pre" in cs.lower() for cs in call_sites)
    has_post = any("post" in cs.lower() for cs in call_sites)

    if DIAG_CALLSITE_MODE == "pre_only":
        if has_pre and not has_post:
            pass_line("Call-site mode validation (pre_only)")
        else:
            fail(f"Expected pre_only but found: pre={has_pre}, post={has_post}")
    elif DIAG_CALLSITE_MODE == "post_only":
        if has_post and not has_pre:
            pass_line("Call-site mode validation (post_only)")
        else:
            fail(f"Expected post_only but found: pre={has_pre}, post={has_post}")
    else:  # "both"
        if has_pre and has_post:
            pass_line("Call-site mode validation (both)")
        else:
            fail(f"Expected both pre and post but found: pre={has_pre}, post={has_post}")

    # -----------------------------------------------------------------
    # 4) SW_TOA accuracy check
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
    # 5) heating_rate_max finite and nonzero
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
    # 6) Stability check via coefficient of variation
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
    # (Expected if cadence includes pre/post at same step and mode="both".)
    # -----------------------------------------------------------------
    step_time = [(to_int(r, "step"), to_float(r, "time")) for r in rows]
    uniq_step_time = len(set(step_time))
    if uniq_step_time == len(step_time):
        pass_line("Unique (step,time) pairs check")
    else:
        print(
            "  Unique (step,time) pairs check [INFO]\n"
            "    Duplicate (step,time) pairs detected. "
            f"This is expected when diag_callsite_mode='{DIAG_CALLSITE_MODE}' with pre+post."
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
