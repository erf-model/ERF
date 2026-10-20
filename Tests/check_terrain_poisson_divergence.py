#!/usr/bin/env python3
"""Check that the anelastic projection converged on one level of a run.

ERF prints, with erf.mg_v >= 1, the largest divergence of the momentum before every
solve (one line per subdomain) and after the solve (one line per level):

  Max/L2 norm of divergence before solve in subdomain 0 at level 1 : 0.44 19.5 and ...
  Max/L2 norm of divergence after  solve at level 1 : 2.4e-08 3.1e-07 and ...

When every face of the region solved on is a wall or a coarse/fine face, the problem is
singular and ERF first subtracts the mean of the divergence (" Subtracting m from rhs");
a projection cannot remove that constant, so a converged solve leaves |m| behind.  For
every solve on the requested level this script takes the residual

    r = | largest divergence after the solve - |m| |        (m = 0 when nothing is subtracted)

and requires r <= max_ratio * (largest divergence before the solve) + atol.  A terrain
Poisson operator that lets flux through a Neumann face of a sloping region cannot converge
there and leaves r of order 0.01 to 1 of the divergence before the solve.

The log gives only the largest values, so r sees the part of the remaining error that adds
to |m|, and on a level with several regions m is the largest of their means; a solve that
fails as the old operator does leaves r three or more orders of magnitude above the limit.

With --max-mean the mean itself is limited too: it is the net inflow through the region's
faces divided by its volume, which is zero when the fluxes through a refined level's
coarse/fine faces add up to the coarse level's (a level that stops inside the domain in z
leaves it nonzero unless the flux through its top and bottom comes from the coarse level).

Exit status: 0 when every solve passes, 1 when one does not, 2 when the log has no solve
on that level (or fewer than --min-solves), or a solve whose divergence before it was
already below --min-before (the check would say nothing).
"""

import argparse
import re
import statistics
import sys

BEFORE = re.compile(r"divergence before solve in subdomain \d+ at level (\d+) : (\S+)")
MEAN = re.compile(r"^\s*Subtracting (\S+) from rhs in subdomain")
AFTER = re.compile(r"divergence after\s+solve at level (\d+) : (\S+)")


def read_solves(path, level):
    """Return (before, mean, after) of every solve on the level: the largest divergence
    before it, the largest mean subtracted from a subdomain, the largest after it."""
    solves = []
    before = None
    mean = 0.0
    current = None  # level of the last 'before solve' line
    with open(path, encoding="utf-8", errors="replace") as log:
        for line in log:
            m = BEFORE.search(line)
            if m:
                current = int(m.group(1))
                if current == level:
                    value = float(m.group(2))
                    before = value if before is None else max(before, value)
                continue
            m = MEAN.search(line)
            if m and current == level:
                mean = max(mean, abs(float(m.group(1))))
                continue
            m = AFTER.search(line)
            if m and int(m.group(1)) == level:
                if before is None:
                    raise ValueError(f"{path}: an 'after solve' line at level {level} "
                                     "without a 'before solve' line")
                solves.append((before, mean, float(m.group(2))))
                before = None
                mean = 0.0
    return solves


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("log", help="ERF run log written with erf.mg_v >= 1")
    parser.add_argument("--level", type=int, required=True, help="AMR level to check")
    parser.add_argument("--max-ratio", type=float, default=1.0e-5,
                        help="largest allowed residual / divergence before a solve")
    parser.add_argument("--atol", type=float, default=2.0e-8,
                        help="residual allowed on top: twice the solvers' default absolute "
                             "tolerance, which multigrid applies to dJ times the divergence")
    parser.add_argument("--min-solves", type=int, default=1,
                        help="fewest solves the log must hold on the level")
    parser.add_argument("--min-before", type=float, default=1.0e-6,
                        help="smallest divergence before a solve for it to count as a check")
    parser.add_argument("--max-mean", type=float, default=None,
                        help="largest allowed |mean| subtracted before a solve (no limit if unset)")
    args = parser.parse_args(argv)

    try:
        solves = read_solves(args.log, args.level)
    except (OSError, ValueError) as err:
        print(f"check_terrain_poisson_divergence: {err}")
        return 2

    if len(solves) < args.min_solves:
        print(f"check_terrain_poisson_divergence: {len(solves)} solves at level {args.level} "
              f"in {args.log}, need at least {args.min_solves}")
        return 2
    small = [b for b, _, _ in solves if b < args.min_before]
    if small:
        print(f"check_terrain_poisson_divergence: a solve at level {args.level} starts from a "
              f"divergence of {min(small):.3e}, below {args.min_before:.1e}; "
              "the check would say nothing")
        return 2

    residuals = [abs(a - m) for _, m, a in solves]
    ratios = [r / b for r, (b, _, _) in zip(residuals, solves)]
    failed = [i for i, (r, (b, _, _)) in enumerate(zip(residuals, solves))
              if r > args.max_ratio * b + args.atol]
    print(f"level {args.level}: {len(solves)} solves; residual / divergence before: "
          f"median {statistics.median(ratios):.3e}, worst {max(ratios):.3e}; "
          f"residual: median {statistics.median(residuals):.3e}, worst {max(residuals):.3e}; "
          f"divergence after: worst {max(a for _, _, a in solves):.3e}, "
          f"mean subtracted: worst {max(m for _, m, _ in solves):.3e}; "
          f"allowed {args.max_ratio:.1e} x before + {args.atol:.1e}")
    if failed:
        i = failed[0]
        b, m, a = solves[i]
        print(f"check_terrain_poisson_divergence: FAIL, {len(failed)} of {len(solves)} solves at "
              f"level {args.level} did not converge; the first (solve {i + 1}) left {residuals[i]:.3e} "
              f"above the mean {m:.3e}, from {b:.3e} before it")
        return 1
    if args.max_mean is not None:
        worst = max(range(len(solves)), key=lambda i: solves[i][1])
        if solves[worst][1] > args.max_mean:
            print(f"check_terrain_poisson_divergence: FAIL, solve {worst + 1} at level {args.level} "
                  f"subtracted a mean of {solves[worst][1]:.3e}, above --max-mean {args.max_mean:.1e}: "
                  "the fluxes through the region's faces do not add up to zero")
            return 1
    print("check_terrain_poisson_divergence: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
