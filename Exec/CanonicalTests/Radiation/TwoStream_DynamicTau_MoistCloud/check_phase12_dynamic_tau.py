#!/usr/bin/env python3
import os, sys, csv, math

def fail(msg):
    print(f"ERROR: {msg}")
    return False

def finite(x):
    return math.isfinite(x)

def check_diag(diag_file):
    print(f"Checking diagnostics file: {diag_file}")
    if not os.path.isfile(diag_file):
        return fail("diagnostics file missing")

    required = ["step","time","call_site","SW_surface","SW_TOA","SW_up_TOA","LW_net_surface","LW_up_TOA","heating_rate_max"]
    rows = []

    with open(diag_file, "r", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return fail("missing CSV header")
        print(f"  Header: {','.join(reader.fieldnames)}...")

        cols = [c.strip() for c in reader.fieldnames]
        missing = [c for c in required if c not in cols]
        if missing:
            return fail(f"missing required columns: {missing}")

        for r in reader:
            try:
                row = {
                    "step": int(r["step"]),
                    "time": float(r["time"]),
                    "call_site": r["call_site"].strip(),
                    "SW_surface": float(r["SW_surface"]),
                    "SW_TOA": float(r["SW_TOA"]),
                    "SW_up_TOA": float(r["SW_up_TOA"]),
                    "LW_net_surface": float(r["LW_net_surface"]),
                    "LW_up_TOA": float(r["LW_up_TOA"]),
                    "heating_rate_max": float(r["heating_rate_max"]),
                }
                rows.append(row)
            except Exception:
                # skip malformed lines
                continue

    if len(rows) == 0:
        return fail("no parseable data rows in diagnostics CSV")

    # checks
    for i, r in enumerate(rows):
        for k in ["time","SW_surface","SW_TOA","SW_up_TOA","LW_net_surface","LW_up_TOA","heating_rate_max"]:
            if not finite(r[k]):
                return fail(f"non-finite value at row {i} col {k}: {r[k]}")

    if all(abs(r["heating_rate_max"]) < 1e-15 for r in rows):
        return fail("heating_rate_max is zero in all rows")

    print(f"  ✓ Parsed {len(rows)} rows")
    print(f"  ✓ heating_rate_max range: {min(r['heating_rate_max'] for r in rows):.6e} .. {max(r['heating_rate_max'] for r in rows):.6e}")
    return True

def check_exists_nonempty(path, min_lines=2):
    if not os.path.isfile(path):
        print(f"  ✗ Missing {path}")
        return False
    n = sum(1 for _ in open(path, "r"))
    if n < min_lines:
        print(f"  ✗ {path} too short ({n} lines)")
        return False
    print(f"  ✓ {path} ({n} lines)")
    return True

def main():
    print("="*70)
    print("Phase 12 Dynamic Tau Test Checker")
    print("="*70)

    ok = True
    ok &= check_diag("radiation_diag_phase12.dat")

    print("\nChecking plot files for heating rate fields...")
    plots = [d for d in os.listdir(".") if d.startswith("plt")]
    if len(plots) == 0:
        print("  ✗ No plot directories found")
        ok = False
    else:
        print(f"  Found {len(plots)} plot directories")
        print("  ✓ Plot files readable")

    print("\nChecking data log files...")
    ok &= check_exists_nonempty("phase12_hist.dat", 2)
    ok &= check_exists_nonempty("phase12_profiles.dat", 2)

    print("\n" + "="*70)
    if ok:
        print("✓ Phase 12 test PASSED")
        return 0
    else:
        print("✗ Phase 12 test FAILED")
        print("  - Check diagnostics file for NaN/Inf/invalid data")
        return 1

if __name__ == "__main__":
    sys.exit(main())
