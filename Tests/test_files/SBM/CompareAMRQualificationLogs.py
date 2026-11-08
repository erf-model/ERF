#!/usr/bin/env python3
"""Compare composite leaf inventories for identical SBM qualification decks."""

import math
import re
import sys


def read_inventory(path):
    values = {}
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            match = re.search(
                r"SBM_AMR_QUAL_INVENTORY component=(\d+) initial=([^ ]+) "
                r"final=([^ ]+).*?abs_error=([^ ]+).*?magnitude_scale=([^ ]+) "
                r"relative_error=([^ ]+) tolerance=([^ ]+) "
                r"leaf_cell_terms=(\d+) operation_count=(\d+) "
                r"reduction_depth=(\d+)",
                line,
            )
            if match:
                component = int(match.group(1))
                values[component] = {
                    "initial": float(match.group(2)),
                    "final": float(match.group(3)),
                    "abs_error": float(match.group(4)),
                    "magnitude_scale": float(match.group(5)),
                    "relative_error": float(match.group(6)),
                    "tolerance": float(match.group(7)),
                    "leaf_cell_terms": int(match.group(8)),
                    "operation_count": int(match.group(9)),
                    "reduction_depth": int(match.group(10)),
                }
    return values


def main():
    if len(sys.argv) != 4:
        raise SystemExit("usage: CompareAMRQualificationLogs.py REF TEST NCOMP")
    reference, candidate, ncomp = sys.argv[1], sys.argv[2], int(sys.argv[3])
    ref = read_inventory(reference)
    got = read_inventory(candidate)
    if set(ref) != set(range(ncomp)) or set(got) != set(range(ncomp)):
        raise SystemExit(
            f"expected exactly {ncomp} inventory components; "
            f"reference={sorted(ref)}, candidate={sorted(got)}"
        )

    eps = sys.float_info.epsilon
    max_inventory_ratio = 0.0
    max_parity_ratio = 0.0
    leaf_cell_terms = None
    for component in range(ncomp):
        reference_row = ref[component]
        candidate_row = got[component]
        if leaf_cell_terms is None:
            leaf_cell_terms = reference_row["leaf_cell_terms"]
        if reference_row["leaf_cell_terms"] != leaf_cell_terms or \
                candidate_row["leaf_cell_terms"] != leaf_cell_terms:
            raise SystemExit("composite leaf-cell count changed between components or runs")

        for label, row in (("reference", reference_row), ("candidate", candidate_row)):
            expected_operations = (
                8 * row["leaf_cell_terms"] + 2 * row["reduction_depth"] + 1
            )
            if row["operation_count"] != expected_operations:
                raise SystemExit(
                    f"component {component} {label} operation count is inconsistent: "
                    f"{row['operation_count']} vs {expected_operations}"
                )
            if row["abs_error"] > row["tolerance"]:
                raise SystemExit(
                    f"component {component} {label} inventory exceeds its analytic bound: "
                    f"error={row['abs_error']:.17g}, bound={row['tolerance']:.17g}"
                )
            ratio = row["abs_error"] / row["tolerance"] if row["tolerance"] > 0 else 0.0
            max_inventory_ratio = max(max_inventory_ratio, ratio)

        parity_operations = max(
            reference_row["operation_count"], candidate_row["operation_count"]
        )
        q = parity_operations * eps
        gamma = q / (1.0 - q)
        for moment in ("initial", "final"):
            ref_value = reference_row[moment]
            got_value = candidate_row[moment]
            bound = gamma * max(
                abs(ref_value) + abs(got_value), sys.float_info.min
            )
            error = abs(ref_value - got_value)
            ratio = error / bound if bound > 0 else 0.0
            max_parity_ratio = max(max_parity_ratio, ratio)
            if error > bound:
                raise SystemExit(
                    f"component {component} {moment} inventory differs: "
                    f"error={error:.17g}, bound={bound:.17g}, gamma={gamma:.3g}"
                )

    print(
        f"component inventories match: leaf_cell_terms={leaf_cell_terms}, "
        f"max_inventory_error_over_bound={max_inventory_ratio:.6g}, "
        f"max_rank_parity_error_over_bound={max_parity_ratio:.6g}"
    )


if __name__ == "__main__":
    main()
