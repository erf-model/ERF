#!/usr/bin/env python3
"""Check that the two-stream prognostic surface state is averaged down.

ERF gives every level its own force-restore surface temperature, and that
temperature is the longwave boundary condition for the column sweep, so the
levels have to agree about the ground they share.  After the finer levels
advance, ERF averages t_sfc and q_sfc down onto the coarse level.

That average-down establishes an exact relation: a coarse cell covered by the
fine level holds the mean of the fine cells above it.  This script asserts that
relation, which makes the tolerance a round-off tolerance rather than a
physical one -- the check separates a correct run (agreement to ~1e-13 K) from
one with no average-down (disagreement at the discretization error, some 1e-4 K
here) by many orders of magnitude, and it does so from the first output rather
than waiting for the levels to drift apart.

Comparing domain means would NOT work: average_down is mean-preserving, so the
coarse and fine means agree whether or not the transfer ever runs.

The surface field in the companion test case varies in x only, so a single
1-d slice from amrex_fextract carries the whole field.  Nothing here assumes
otherwise, but a case with y-structure would need the full field instead.
"""

import argparse
import subprocess
import sys


def parse_slice(path):
    """Read an amrex_fextract slice file into (coords, values)."""
    coords, values = [], []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}: expected 'x value', got {line!r}")
            coords.append(float(fields[0]))
            values.append(float(fields[1]))
    if not values:
        raise ValueError(f"{path}: no data rows")
    return coords, values


def extract(fextract, plotfile, variable, level, out_path, direction=0):
    """Run amrex_fextract for a single level and return its slice."""
    cmd = [fextract, '-d', str(direction), '-v', variable,
           '-c', str(level), '-f', str(level),
           '-e', '-s', out_path, plotfile]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"amrex_fextract failed for level {level} (exit {proc.returncode})\n"
            f"  command: {' '.join(cmd)}\n{proc.stdout}{proc.stderr}")
    return parse_slice(out_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--plotfile', required=True,
                        help='2-d plotfile holding the surface state')
    parser.add_argument('--fextract', required=True,
                        help='path to the amrex_fextract executable')
    parser.add_argument('--var', default='seb_t_sfc',
                        help='variable to check (default: seb_t_sfc)')
    parser.add_argument('--coarse-level', type=int, default=0)
    parser.add_argument('--fine-level', type=int, default=1)
    parser.add_argument('--ref-ratio', type=int, default=2,
                        help='refinement ratio in the slice direction')
    parser.add_argument('--tol', type=float, default=1.0e-8,
                        help='absolute tolerance; a round-off bound, not a '
                             'physical one (default: 1e-8)')
    parser.add_argument('--scratch-prefix', default='seb_parity_slice',
                        help='prefix for the slice files this writes')
    args = parser.parse_args()

    if args.ref_ratio < 2:
        parser.error('--ref-ratio must be at least 2')
    if args.fine_level <= args.coarse_level:
        parser.error('--fine-level must be above --coarse-level')
    if args.tol <= 0.0:
        parser.error('--tol must be positive')

    try:
        _, coarse = extract(args.fextract, args.plotfile, args.var,
                            args.coarse_level, f'{args.scratch_prefix}_c.txt')
        _, fine = extract(args.fextract, args.plotfile, args.var,
                          args.fine_level, f'{args.scratch_prefix}_f.txt')
    except (OSError, RuntimeError, ValueError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2

    ratio = args.ref_ratio ** (args.fine_level - args.coarse_level)
    if len(fine) != ratio * len(coarse):
        print(f'ERROR: level {args.fine_level} has {len(fine)} points but level '
              f'{args.coarse_level} has {len(coarse)}; expected a factor of '
              f'{ratio}. The fine level must cover the whole domain for this '
              f'check to be meaningful.', file=sys.stderr)
        return 2

    # A field that never leaves its initial value would satisfy the invariant
    # trivially, so refuse to report success on one.
    spread = max(fine) - min(fine)
    if spread <= args.tol:
        print(f'ERROR: {args.var} varies by only {spread:.3e} across level '
              f'{args.fine_level}, at or below the tolerance {args.tol:.3e}. '
              f'The field is uniform, so the check would pass on any transfer '
              f'-- and prove nothing. Fix the case, not the tolerance.',
              file=sys.stderr)
        return 2

    worst = 0.0
    worst_i = 0
    for i, c_val in enumerate(coarse):
        block = fine[i * ratio:(i + 1) * ratio]
        expected = sum(block) / float(ratio)
        err = abs(c_val - expected)
        if err > worst:
            worst, worst_i = err, i

    print(f'{args.var}: {len(coarse)} coarse cells vs level {args.fine_level}, '
          f'ref_ratio {ratio}')
    print(f'  field spread on the fine level : {spread:.6e}')
    print(f'  worst |coarse - mean(fine)|    : {worst:.6e} at i = {worst_i}')
    print(f'  tolerance                      : {args.tol:.6e}')

    if worst > args.tol:
        print(f'FAIL: level {args.coarse_level} does not hold the average of '
              f'level {args.fine_level}. The prognostic surface state was not '
              f'averaged down, so the levels disagree about the ground they '
              f'share.', file=sys.stderr)
        return 1

    print('PASS: the coarse surface state is the average of the fine one.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
