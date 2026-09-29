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

With --created-from the script checks the other transfer instead: the one that
gives a fine level created mid-run its starting surface.  ERF interpolates it
from the parent (ERF::fill_seb_from_coarse), so the new level carries on from
the surface the coarse level had reached.  Without that it would start from the
erf.rad_t_sfc scalar, and the average-down would then drag the coarse surface
under it back there too.  See check_created_from for the assertion.

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


def finest_level(plotfile):
    """Return the finest level recorded in a native plotfile's Header."""
    with open(f'{plotfile}/Header') as handle:
        lines = [line.strip() for line in handle]
    # version, ncomp, ncomp variable names, ndim, time, finest level
    ncomp = int(lines[1])
    return int(lines[2 + ncomp + 2])


def check_created_from(args, fine, ratio):
    """Assert that a fine level created mid-run started from its parent's surface.

    Takes three 2-d plotfiles written while only the coarse level existed: the
    initial one, and the last two before the fine level was created (PREV, BEFORE).
    --plotfile is the first one written after the creation, one step past BEFORE.

    Between creation and that output the new level advances one step, so it holds
    its starting surface plus one step's change.  The coarse surface drifts steadily,
    so one step past BEFORE the parent's surface is BEFORE + (BEFORE - PREV), and a
    level that started from the parent must average to that over each coarse cell.
    What separates the two transfers is how far BEFORE sits from the initial surface:
    a level built from the scalar lands that far away instead.  So the tolerance is a
    fraction of one step's change, and the check refuses to run unless the drift at
    creation is several tolerances wide.
    """
    names = args.created_from.split(',')
    if len(names) != 3:
        print(f'ERROR: --created-from takes INITIAL,PREV,BEFORE; got {args.created_from!r}',
              file=sys.stderr)
        return 2
    # The whole point is a level that did NOT exist yet. If it was already there, the
    # check below is comparing a level with its own history, not testing its creation.
    for name in names:
        if finest_level(name) != args.coarse_level:
            print(f'ERROR: {name} already holds levels above {args.coarse_level}. '
                  f'The fine level must be created after these outputs, or its '
                  f'creation is not what is being tested.', file=sys.stderr)
            return 2
    try:
        initial, prev, before = (
            extract(args.fextract, name, args.var, args.coarse_level,
                    f'{args.scratch_prefix}_created{n}.txt')[1]
            for n, name in enumerate(names))
    except (OSError, RuntimeError, ValueError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2

    step = max(abs(b - p) for b, p in zip(before, prev))
    tol = args.created_tol_frac * step
    # A level built from the scalar would start at the initial surface, so it misses
    # by the drift BEFORE has accumulated since then.
    drift = min(abs(b - i) for b, i in zip(before, initial))
    print(f'{args.var}: level {args.fine_level} created after {names[2]}')
    print(f'  largest one-step change on level {args.coarse_level}: {step:.6e}')
    print(f'  smallest drift from {names[0]}           : {drift:.6e}')
    print(f'  tolerance ({args.created_tol_frac} of one step)       : {tol:.6e}')
    if step <= 0.0 or drift <= 4.0 * tol:
        print(f'ERROR: the coarse surface had drifted only {drift:.3e} from its initial '
              f'value when the level was created, against a tolerance of {tol:.3e}. A '
              f'level started from the scalar would pass too. Create the level later, '
              f'or drive the surface harder -- not a looser tolerance.', file=sys.stderr)
        return 2

    worst = 0.0
    worst_i = 0
    for i, (b_val, p_val) in enumerate(zip(before, prev)):
        expected = 2.0 * b_val - p_val
        got = sum(fine[i * ratio:(i + 1) * ratio]) / float(ratio)
        err = abs(got - expected)
        if err > worst:
            worst, worst_i = err, i
    print(f'  worst |mean(fine) - extrapolated parent|: {worst:.6e} at i = {worst_i}')

    if worst > tol:
        print(f'FAIL: level {args.fine_level} did not start from the surface level '
              f'{args.coarse_level} had reached when it was created. The new level was '
              f'filled with something else (the erf.rad_t_sfc scalar, if the transfer '
              f'is missing).', file=sys.stderr)
        return 1
    print('PASS: the new level started from its parent\'s surface.')
    return 0


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
    parser.add_argument('--evolved-from', type=float, default=None,
                        help='instead of comparing levels, assert every cell of the COARSE '
                             'level has moved away from this value (the initial '
                             'erf.rad_t_sfc). Catches a coarse surface pinned at its '
                             'level-creation value by an average-down off a level that '
                             'never advanced its own copy.')
    parser.add_argument('--created-from', default=None, metavar='INITIAL,PREV,BEFORE',
                        help='instead of comparing levels, assert the fine level in '
                             '--plotfile, created one step after BEFORE, started from the '
                             'coarse surface rather than the erf.rad_t_sfc scalar. The '
                             'three 2-d plotfiles hold the coarse level only.')
    parser.add_argument('--created-tol-frac', type=float, default=0.25,
                        help='tolerance for --created-from, as a fraction of the '
                             'largest one-step change of the coarse surface (default: 0.25)')
    args = parser.parse_args()

    if args.ref_ratio < 2:
        parser.error('--ref-ratio must be at least 2')
    if args.fine_level <= args.coarse_level:
        parser.error('--fine-level must be above --coarse-level')
    if args.tol <= 0.0:
        parser.error('--tol must be positive')
    if args.created_from is not None and args.evolved_from is not None:
        parser.error('--created-from and --evolved-from are separate checks')
    if not 0.0 < args.created_tol_frac < 1.0:
        parser.error('--created-tol-frac must lie strictly between 0 and 1')

    try:
        _, coarse = extract(args.fextract, args.plotfile, args.var,
                            args.coarse_level, f'{args.scratch_prefix}_c.txt')
        _, fine = extract(args.fextract, args.plotfile, args.var,
                          args.fine_level, f'{args.scratch_prefix}_f.txt')
    except (OSError, RuntimeError, ValueError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2

    # A shallow nest does not sweep, so its surface state stays frozen at what
    # fill_seb_from_coarse wrote and there is no fine solution to compare against. What
    # must hold instead is that the coarse level kept evolving underneath it.
    if args.evolved_from is not None:
        worst_cell = min(coarse, key=lambda v: abs(v - args.evolved_from))
        gap = abs(worst_cell - args.evolved_from)
        print(f'{args.var}: {len(coarse)} cells on level {args.coarse_level}')
        print(f'  closest cell to {args.evolved_from}: {worst_cell!r} (gap {gap:.6e})')
        print(f'  tolerance                         : {args.tol:.6e}')
        if gap <= args.tol:
            print(f'FAIL: a cell of level {args.coarse_level} is still at '
                  f'{args.evolved_from}, its level-creation value. The coarse surface was '
                  f'overwritten by a level that never advanced its own copy.',
                  file=sys.stderr)
            return 1
        print('PASS: the coarse surface state evolved everywhere.')
        return 0

    ratio = args.ref_ratio ** (args.fine_level - args.coarse_level)
    if len(fine) != ratio * len(coarse):
        print(f'ERROR: level {args.fine_level} has {len(fine)} points but level '
              f'{args.coarse_level} has {len(coarse)}; expected a factor of '
              f'{ratio}. The fine level must cover the whole domain for this '
              f'check to be meaningful.', file=sys.stderr)
        return 2

    if args.created_from is not None:
        return check_created_from(args, fine, ratio)

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
