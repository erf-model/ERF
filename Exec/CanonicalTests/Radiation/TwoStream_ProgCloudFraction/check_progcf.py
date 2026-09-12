#!/usr/bin/env python3
"""
Regression Test Checker: Prognostic Cloud Fraction for TwoStream Radiation

Validates:
1. The radiation diagnostics CSV exists and carries the expected columns
2. Every radiative flux and heating rate is finite
3. Shortwave fluxes are physically ordered (0 <= SW_surface <= SW_TOA,
   0 <= SW_up_TOA <= SW_TOA)
4. The column heats: heating_rate_max is non-zero
"""

import csv
import math
import os
import sys

DIAG_FILE = "radiation_progcf_diag.dat"

# Columns written for every run. The SEB columns that follow are NaN by design
# unless the surface energy balance is enabled, so they are not checked here.
REQUIRED = [
    "step", "time", "call_site", "SW_surface", "SW_TOA", "SW_up_TOA",
    "LW_net_surface", "LW_up_TOA", "heating_rate_max",
]
NUMERIC = REQUIRED[3:]


def fail(msg):
    print(f"ERROR: {msg}")
    return False


def read_rows(path):
    """Return (rows, error). A missing or empty file is an error, not a pass."""
    if not os.path.isfile(path):
        return None, f"diagnostics file not found: {path}"
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return None, f"no CSV header in {path}"
        cols = [c.strip() for c in reader.fieldnames]
        missing = [c for c in REQUIRED if c not in cols]
        if missing:
            return None, f"missing required columns in {path}: {missing}"
        rows = [r for r in reader if r.get("step")]
    if not rows:
        return None, f"no data rows in {path}"
    return rows, None


def check_diag(path):
    rows, err = read_rows(path)
    if err:
        return fail(err)

    heating = []
    for i, r in enumerate(rows):
        vals = {}
        for c in NUMERIC:
            try:
                vals[c] = float(r[c])
            except (TypeError, ValueError):
                return fail(f"row {i}: column {c} is not a number: {r[c]!r}")
            if not math.isfinite(vals[c]):
                return fail(f"row {i}: column {c} is not finite: {vals[c]}")

        if vals["SW_TOA"] < 0.0:
            return fail(f"row {i}: negative SW_TOA {vals['SW_TOA']}")
        if not (0.0 <= vals["SW_surface"] <= vals["SW_TOA"] + 1.0e-6):
            return fail(
                f"row {i}: SW_surface {vals['SW_surface']} outside "
                f"[0, SW_TOA={vals['SW_TOA']}]"
            )
        if not (0.0 <= vals["SW_up_TOA"] <= vals["SW_TOA"] + 1.0e-6):
            return fail(
                f"row {i}: SW_up_TOA {vals['SW_up_TOA']} outside "
                f"[0, SW_TOA={vals['SW_TOA']}]"
            )
        heating.append(vals["heating_rate_max"])

    if all(abs(h) < 1.0e-15 for h in heating):
        return fail("heating_rate_max is zero in every row; radiation did not heat the column")

    print(f"  Parsed {len(rows)} rows from {os.path.basename(path)}")
    print(f"  heating_rate_max range: {min(heating):.6e} .. {max(heating):.6e} K/s")
    print(f"  SW_surface range: {min(float(r['SW_surface']) for r in rows):.3f} .. "
          f"{max(float(r['SW_surface']) for r in rows):.3f} W/m^2")
    return True


def check_plotfiles():
    plots = sorted(d for d in os.listdir(".") if d.startswith("plt") and os.path.isdir(d))
    if not plots:
        return fail("no plotfile directories were written")
    for d in plots:
        if not os.path.isfile(os.path.join(d, "Header")):
            return fail(f"plotfile {d} has no Header")
    print(f"  {len(plots)} plotfiles written, all with a Header")
    return True


def main():
    print("=" * 70)
    print("Regression Test: Prognostic Cloud Fraction for TwoStream")
    print("=" * 70)

    ok = check_diag(DIAG_FILE)
    ok = check_plotfiles() and ok

    print("=" * 70)
    if ok:
        print("RESULT: PASS")
        return 0
    print("RESULT: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
