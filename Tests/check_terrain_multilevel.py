#!/usr/bin/env python3
"""Check the projections of one level in an ERF log written with erf.mg_v = 1.

For every solve on the requested level the log holds the divergence before the
solve, the compatibility constant subtracted from the right-hand side when the
problem is singular (every face Neumann or periodic, as on a refined level), and the
divergence after the velocity correction.  A converged solve removes everything but
that constant, so the max norm of the divergence after the solve must equal the
constant's magnitude to within tol times the divergence before the solve.  The exit
code is the verdict; a table of every solve is printed.
"""
import argparse
import re
import sys

BEFORE = re.compile(r"divergence before solve in subdomain \d+ at level (\d+) : ([-+0-9.eE]+) ")
SUBTRACT = re.compile(r"Subtracting ([-+0-9.eE]+) from rhs in subdomain")
AFTER = re.compile(r"divergence after  solve at level (\d+) : ([-+0-9.eE]+) ")


def parse(path, level):
    solves = []
    before = None
    offset = 0.0
    for line in open(path, errors="replace"):
        m = BEFORE.search(line)
        if m:
            before = (int(m.group(1)), float(m.group(2)))
            offset = 0.0
            continue
        m = SUBTRACT.search(line)
        if m and before is not None:
            offset = float(m.group(1))
            continue
        m = AFTER.search(line)
        if m and before is not None and int(m.group(1)) == before[0]:
            if before[0] == level:
                solves.append((before[1], offset, float(m.group(2))))
            before = None
    return solves


def verdict(solves, level, min_solves, tol):
    """Return (exit code, table lines) for the parsed solves of one level."""
    lines = []
    if len(solves) < min_solves:
        lines.append(f"found {len(solves)} projections on level {level}, fewer than {min_solves}")
        return 1, lines
    bad = 0
    lines.append(f"{'before':>12} {'offset':>12} {'after':>12} {'|after|-|offset|':>18} {'allowed':>12}")
    for before, offset, after in solves:
        left = abs(abs(after) - abs(offset))
        allowed = tol * before + 1.0e-14
        flag = "" if left <= allowed else "  FAIL"
        if flag:
            bad += 1
        lines.append(f"{before:12.4e} {offset:12.4e} {after:12.4e} {left:18.4e} {allowed:12.4e}{flag}")
    lines.append(f"{len(solves)} projections on level {level}, {bad} not converged")
    return (1 if bad else 0), lines


def self_test():
    """The checker's own logic on synthetic logs: converged passes, 1e-3 off fails, missing level fails."""
    import os
    import tempfile
    head = ("Max/L2 norm of divergence before solve in subdomain 0 at level 1 : 6.0e-06 1.0e-05 and volume-weighted sum 3\n"
            " Subtracting 9.1765e-08 from rhs in subdomain 0\n"
            "Solving the terrain Poisson equation with MLTerrainPoisson multigrid\n")
    good = head + "Max/L2 norm of divergence after  solve at level 1 : 9.1767e-08 1.4e-05 and volume-weighted sum 348\n"
    bad = head + "Max/L2 norm of divergence after  solve at level 1 : 9.7e-08 1.4e-05 and volume-weighted sum 348\n"
    other = good.replace("level 1", "level 0")
    failures = 0
    for name, text, expect in (("converged", good, 0), ("off by 1e-3 of before", bad, 1), ("no level-1 report", other, 1)):
        with tempfile.NamedTemporaryFile("w", suffix=".log", delete=False) as f:
            f.write(text)
            path = f.name
        code, lines = verdict(parse(path, 1), 1, 1, 1.0e-6)
        os.unlink(path)
        ok = (code == expect)
        print(f"self-test {name}: exit {code}, expected {expect}: {'ok' if ok else 'WRONG'}")
        if not ok:
            failures += 1
    return 1 if failures else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", nargs="?")
    ap.add_argument("--level", type=int)
    ap.add_argument("--min-solves", type=int, default=1)
    ap.add_argument("--tol", type=float, default=1.0e-6,
                    help="allowed |after - |offset|| as a fraction of the divergence before the solve")
    ap.add_argument("--self-test", action="store_true", help="check the verdict logic on synthetic logs")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if args.log is None or args.level is None:
        ap.error("a log and --level are required")

    solves = parse(args.log, args.level)
    code, lines = verdict(solves, args.level, args.min_solves, args.tol)
    print("\n".join(lines))
    return code


if __name__ == "__main__":
    sys.exit(main())
