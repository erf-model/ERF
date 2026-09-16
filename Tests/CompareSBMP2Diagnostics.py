#!/usr/bin/env python3
"""Compare scalar production SBM qualification diagnostics across MPI layouts."""

import math
import sys


def read(path):
    values = {}
    with open(path, encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key] = value
    return values


def numeric_keys(values):
    fixed = {
        "accepted_face_transfer_l1",
        "accepted_face_transfer_max",
        "accepted_bulk_transfer_l1",
        "accepted_bulk_transfer_max",
        "minimum_accepted_limiter",
        "composite_initial_total",
        "composite_final_total",
        "compact_initial_total",
        "compact_final_total",
        "compact_final_qc",
        "compact_final_qr",
        "spectral_cloud_initial",
        "spectral_cloud_final",
        "spectral_rain_initial",
        "spectral_rain_final",
        "accepted_bulk_transfer_sum_qc",
        "accepted_bulk_transfer_sum_qr",
    }
    fixed.update(key for key in values if key.startswith("composite_initial_comp_"))
    fixed.update(key for key in values if key.startswith("composite_final_comp_"))
    fixed.update(key for key in values if key.startswith("accepted_face_transfer_sum_comp_"))
    return sorted(fixed)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: CompareSBMP2Diagnostics.py one_rank two_rank")
    one = read(sys.argv[1])
    two = read(sys.argv[2])
    missing = sorted(set(numeric_keys(one)) - set(two))
    if missing:
        raise SystemExit("two-rank diagnostic is missing: " + ", ".join(missing))
    for key in numeric_keys(one):
        try:
            first = float(one[key])
            second = float(two[key])
        except ValueError as error:
            raise SystemExit(f"non-numeric comparison field {key}: {error}") from error
        if not math.isclose(first, second, rel_tol=5.0e-13, abs_tol=5.0e-18):
            raise SystemExit(f"MPI decomposition mismatch for {key}: {first} vs {second}")
    print("1-rank/2-rank accepted transfer, state, projection, and total diagnostics agree")


if __name__ == "__main__":
    main()
