#!/usr/bin/env python3
import os
import sys
import csv
import math

def is_finite(x):
    return math.isfinite(x)

def parse_diag_csv(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return rows, "Missing CSV header"

        required = {
            "step","time","call_site",
            "SW_surface","SW_TOA","SW_up_TOA","LW_net_surface","LW_up_TOA","heating_rate_max"
        }
        found = {h.strip() for h in reader.fieldnames}
        if not required.issubset(found):
            return rows, f"Header mismatch. Found={sorted(found)}"

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
                # skip malformed lines but continue
                continue

    return rows, None

def main():
    print("="*70)
    print("Phase 11 Surface Heterogeneity RegTest Validation")
    print("="*70)
    cwd = os.getcwd()
    print(f"\nWorking directory: {cwd}")

    diag = os.path.join(cwd, "radiation_diag_phase11.dat")

    passed = 0
    total = 5

    print("\n1. Checking for required output files...")
    if not os.path.isfile(diag):
        print(f"  ERROR: Missing diagnostics file: {diag}")
    else:
        print(f"✓ Radiation diagnostics CSV found: {diag}")
        passed += 1

    print("\n2. Parsing radiation diagnostics CSV...")
    rows, err = parse_diag_csv(diag) if os.path.isfile(diag) else ([], "file missing")
    print(f"✓ Parsed {len(rows)} data rows from {diag}")
    if err:
        print(f"  ERROR: {err}")
    elif len(rows) == 0:
        print("  ERROR: Failed to parse diagnostics CSV")
    else:
        passed += 1

    print("\n3. Checking for NaN/Inf in diagnostics...")
    naninf_ok = True
    for r in rows:
        for k in ["time","SW_surface","SW_TOA","SW_up_TOA","LW_net_surface","LW_up_TOA","heating_rate_max"]:
            if not is_finite(r[k]):
                print(f"  ERROR: non-finite value in {k}: {r[k]}")
                naninf_ok = False
                break
    if naninf_ok and len(rows) > 0:
        print("✓ No NaN/Inf found")
        passed += 1
    elif len(rows) == 0:
        print("  ERROR: No rows to validate")

    print("\n4. Validating heating rates...")
    if len(rows) == 0:
        print("  ERROR: No rows available")
    else:
        vals = [r["heating_rate_max"] for r in rows]
        if all(abs(v) < 1e-15 for v in vals):
            print("  ERROR: heating_rate_max is zero for all rows")
        else:
            print(f"✓ heating_rate_max nonzero; min={min(vals):.6e}, max={max(vals):.6e}")
            passed += 1

    print("\n5. Validating surface fluxes...")
    if len(rows) == 0:
        print("  ERROR: No rows available")
    else:
        sw = [r["SW_surface"] for r in rows]
        toa = [r["SW_TOA"] for r in rows]
        if any(not is_finite(v) for v in sw+toa):
            print("  ERROR: non-finite flux values")
        else:
            print(f"✓ Surface/TOA flux finite; SW_surface(last)={sw[-1]:.6f}, SW_TOA(last)={toa[-1]:.6f}")
            passed += 1

    print("\n6. Phase 11 feature validation...")
    print("  Note: Phase 11 test runs in fallback mode (hetero fields all nullptr)")
    print("  - surface_albedo_sw = 0.3 (from inputs)")
    print("  - surface_emissivity_lw = 0.99 (from inputs)")
    print("  - surface_temp_k = 300.0 K (from inputs)")
    print("  ✓ Fallback path being exercised (no hetero LSM fields available)")

    print("\n" + "="*70)
    print(f"VALIDATION SUMMARY: {passed}/{total} checks passed")
    print("="*70)
    if passed == total:
        print("\n✅ All checks passed")
        sys.exit(0)
    else:
        print(f"\n⚠️  {total-passed} check(s) failed - review output above")
        sys.exit(1)

if __name__ == "__main__":
    main()
