#!/usr/bin/env python3
"""Check the lines of LineSampling_LatLon.i: each sample's theta names its column."""

import re
import sys


def theta(i, j):
    return 300.0 + 0.1 * i + 0.001 * j


EXPECTED = {
    "prof": [theta(10, 20)] * 8,
    "ew":   [theta(i, 5) for i in range(4, 26)],
    "ns":   [theta(7, j) for j in range(3, 29)],
}


def main():
    log = open(sys.argv[1]).read()
    for name, expected in EXPECTED.items():
        rows = [r.split() for r in open(f"{name}.theta") if r.strip()]
        if not rows:
            raise SystemExit(f"{name}.theta is empty")
        for row in rows:
            values = [float(v) for v in row[1:]]
            if len(values) != len(expected):
                raise SystemExit(f"{name}: {len(values)} samples at t = {row[0]}, expected {len(expected)}")
            for n, (got, want) in enumerate(zip(values, expected)):
                if abs(got - want) > 1.0e-4:
                    raise SystemExit(f"{name}: sample {n} at t = {row[0]} has theta {got}, expected {want}")
        print(f"{name}: {len(rows)} time(s) of {len(expected)} samples match")
    # The hi ends are 0.4 and 0.7 cells of 1000 m across their lines
    for line, snap in (("Line 1 of", 400.0), ("Line 2 of", 700.0)):
        match = [re.search(r"through the lo end, ([0-9.eE+-]+) m from the hi end", l)
                 for l in log.splitlines() if l.startswith(line)]
        match = [m for m in match if m]
        if len(match) != 1 or abs(float(match[0].group(1)) - snap) > 1.0:
            raise SystemExit(f"the log does not report {line[:-3]} {snap} m from its hi end")
        print(f"{line[:-3]}: hi end {float(match[0].group(1)):.1f} m across the line")

    # erf.sample_line_latlon: the column at the centre of cell (10,20), and its log written
    match = re.search(r"erf\.sample_line_latlon line 0 at .* \(x=([0-9.eE+-]+) y=([0-9.eE+-]+)\)", log)
    if not match or abs(float(match.group(1)) - 10500.0) > 1.0 or abs(float(match.group(2)) - 20500.0) > 1.0:
        raise SystemExit("erf.sample_line_latlon was not placed at the centre of cell (10,20): " +
                         (match.group(0) if match else "no position printed"))
    rows = [r for r in open("col_ll") if r.strip() and not r.startswith("#")]
    if not rows:
        raise SystemExit("the log col_ll of erf.sample_line_latlon is empty")
    print(f"col_ll: placed at x={match.group(1)} y={match.group(2)}, {len(rows)} row(s)")


if __name__ == "__main__":
    main()
