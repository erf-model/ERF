#!/usr/bin/env python3
"""Require fields of an AMReX plotfile to lie within physical bounds.

For each --range NAME:LO:HI the field must exist in the plotfile, be finite
everywhere, and have its minimum and maximum inside [LO, HI]. amrex_fextrema
supplies the extrema, so nothing here depends on the plotfile layout.

Written for the Noah-MP smoke test, where the fields are the land model's own
outputs and the bounds are physical rather than tuned: absorbed shortwave cannot
be negative or exceed the solar constant, a land surface is not below 150 K, and
so on. A bound that has to be tuned to make a run pass is the wrong bound.

The field-fill value -999 (MissingPolicy::FillMinus999WhenUnavailable in the 2D
plotfile catalog) is reported as "not provided" rather than as out of range, since
it means the source never wrote the field, which is a different defect.

Exit status: 0 all fields in range, 1 a field out of range or non-finite,
2 the check could not be made (missing plotfile or field, fextrema failure,
malformed --range).
"""

import argparse
import math
import os
import subprocess
import sys

FILL = -999.0


def parse_range(text):
    parts = text.split(':')
    if len(parts) != 3:
        raise ValueError(f"--range must be NAME:LO:HI, got {text!r}")
    name, lo, hi = parts[0], float(parts[1]), float(parts[2])
    if not name:
        raise ValueError(f"--range has an empty field name: {text!r}")
    if lo > hi:
        raise ValueError(f"--range {text!r} has LO > HI")
    return name, lo, hi


def read_extrema(fextrema, plotfile):
    proc = subprocess.run([fextrema, plotfile], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"amrex_fextrema failed on {plotfile} "
                           f"(exit {proc.returncode})\n{proc.stdout}{proc.stderr}")
    extrema = {}
    in_table = False
    for line in proc.stdout.splitlines():
        if 'minimum value' in line and 'maximum value' in line:
            in_table = True
            continue
        if not in_table:
            continue
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            extrema[fields[0]] = (float(fields[1]), float(fields[2]))
        except ValueError:
            continue
    if not extrema:
        raise RuntimeError(f"no variables parsed from amrex_fextrema output for {plotfile}:\n"
                           f"{proc.stdout}")
    return extrema


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--plotfile', required=True)
    parser.add_argument('--fextrema', required=True, help='path to amrex_fextrema')
    parser.add_argument('--range', action='append', default=[], dest='ranges',
                        metavar='NAME:LO:HI', help='repeatable')
    args = parser.parse_args()

    if not args.ranges:
        parser.error('give at least one --range')
    try:
        ranges = [parse_range(r) for r in args.ranges]
    except ValueError as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2

    # The exit status of the run that wrote this cannot be trusted on its own (a Fortran
    # STOP inside a model exits 0), so the plotfile the run should have reached is part
    # of the check.
    if not os.path.isfile(os.path.join(args.plotfile, 'Header')):
        print(f'ERROR: {args.plotfile} was not written -- the run did not reach the step '
              f'that writes it.', file=sys.stderr)
        return 2
    try:
        extrema = read_extrema(args.fextrema, args.plotfile)
    except (OSError, RuntimeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2

    failed = False
    for name, lo, hi in ranges:
        if name not in extrema:
            print(f'ERROR: {name} is not in {args.plotfile} '
                  f'(have: {", ".join(sorted(extrema))})', file=sys.stderr)
            return 2
        vmin, vmax = extrema[name]
        if vmin == FILL and vmax == FILL:
            print(f'FAIL: {name} is -999 everywhere -- not provided by its source')
            failed = True
            continue
        if not (math.isfinite(vmin) and math.isfinite(vmax)):
            print(f'FAIL: {name} is not finite (min {vmin}, max {vmax})')
            failed = True
            continue
        ok = lo <= vmin and vmax <= hi
        print(f'{"ok  " if ok else "FAIL"} {name:<22} min {vmin:14.6g}  max {vmax:14.6g}  '
              f'bounds [{lo:g}, {hi:g}]')
        failed |= not ok

    if failed:
        print('FAIL: at least one field is outside its physical bounds', file=sys.stderr)
        return 1
    print('PASS: every field is within its physical bounds.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
