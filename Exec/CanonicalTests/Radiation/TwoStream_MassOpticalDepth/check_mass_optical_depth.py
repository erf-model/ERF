#!/usr/bin/env python3
"""
Two-Stream Radiation Validation Script
Mass-based optical depth model: resolution independence

inputs_coarse and inputs_fine run the same moist column on 32 and 64 layers
with tau_model = mass. Because every layer's optical depth follows its mass
path, the column optical depth, and hence the surface and top-of-atmosphere
fluxes, must agree between the two grids. (With the per-layer model the
64-layer column would have twice the optical depth of the 32-layer one.)

Checks, on the last diagnostics row of each run:
1. SW_surface, SW_up_TOA and LW_up_TOA agree within 1% (relative)
2. LW_net_surface agrees within 2 W/m^2 (its magnitude is set by the small
   temperature difference between the surface and the first air layer, whose
   height differs between the grids)
3. all values finite, SW_up_TOA > 0 (Rayleigh + surface reflection)
"""
import csv
import math
import os
import sys


def read_last_row(filename):
    if not os.path.exists(filename):
        print(f"ERROR: Diagnostic file {filename} not found")
        return None
    with open(filename, 'r') as f:
        rows = [r for r in csv.DictReader(f) if any((v or '').strip() for v in r.values())]
    if not rows:
        print(f"ERROR: No data in {filename}")
        return None
    return {k.strip(): v.strip() for k, v in rows[-1].items()}


def main():
    coarse = read_last_row("radiation_mass_coarse_diag.dat")
    fine = read_last_row("radiation_mass_fine_diag.dat")
    if coarse is None or fine is None:
        return False

    print("=" * 70)
    print("Two-Stream Radiation: mass-based optical depth, 32 vs 64 layers")
    print("=" * 70)
    errors = []
    for key in ["SW_surface", "SW_up_TOA", "LW_up_TOA", "LW_net_surface"]:
        c = float(coarse[key])
        f = float(fine[key])
        if not (math.isfinite(c) and math.isfinite(f)):
            errors.append(f"{key} not finite: coarse={c} fine={f}")
            print(f"  {key}: coarse={c} fine={f} [FAIL - not finite]")
            continue
        if key == "LW_net_surface":
            diff = abs(c - f)
            ok = diff <= 2.0
            print(f"  {key}: coarse={c:.4f} fine={f:.4f} |diff|={diff:.4f} W/m^2 {'[PASS]' if ok else '[FAIL]'}")
        else:
            rel = abs(c - f) / max(abs(c), 1.0e-12)
            ok = rel <= 1.0e-2
            print(f"  {key}: coarse={c:.4f} fine={f:.4f} rel diff={rel:.3e} {'[PASS]' if ok else '[FAIL]'}")
        if not ok:
            errors.append(f"{key} differs between the grids")

    sw_up = float(coarse["SW_up_TOA"])
    if not sw_up > 0.0:
        errors.append("SW_up_TOA is not positive")
        print("  SW_up_TOA > 0 [FAIL]")
    else:
        print("  SW_up_TOA > 0 [PASS]")

    print("=" * 70)
    if errors:
        print("TEST FAILED")
        for e in errors:
            print("  - " + e)
        return False
    print("TEST PASSED - fluxes are independent of the vertical resolution with tau_model = mass")
    return True


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
