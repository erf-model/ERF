#!/usr/bin/env python3
"""Check the erf.sample_line logs against erf.do_line_sampling text output of the same columns.

Usage: check_sample_line.py LOG:PREFIX [LOG:PREFIX ...]

Each LOG is an erf.sample_line log and PREFIX the erf.sample_line_name of a do_line_sampling
line along the same column, sampled with the variables erf.sample_line writes.  Every row of
LOG must hold the time and then, in order, each conserved component, the three velocities and
the six stresses over the column, and equal the do_line_sampling values at the same time.
"""
import glob
import sys

VARS_AFTER_CONS = ["x_velocity", "y_velocity", "z_velocity",
                   "sgs_tau11", "sgs_tau12", "sgs_tau13", "sgs_tau22", "sgs_tau23", "sgs_tau33"]
# The log and the text files write 6 significant digits
RTOL, ATOL = 2.0e-5, 1.0e-12


def read_rows(path):
    with open(path) as f:
        return [[float(x) for x in line.split()] for line in f if line.strip()]


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    failed = 0
    for pair in argv[1:]:
        log, prefix = pair.split(":")
        ncons = len(glob.glob(prefix + ".cons_*"))
        if ncons == 0:
            print("%s: no %s.cons_* files" % (log, prefix))
            failed += 1
            continue
        varnames = ["cons_%d" % n for n in range(ncons)] + VARS_AFTER_CONS
        ref = {name: {round(r[0], 9): r[1:] for r in read_rows(prefix + "." + name)}
               for name in varnames}
        nz = len(next(iter(ref[varnames[0]].values())))
        nvals = 1 + nz * len(varnames)
        rows = read_rows(log)
        nchecked = 0
        for row in rows:
            if len(row) != nvals:
                print("%s: a row has %d values, expected %d" % (log, len(row), nvals))
                failed += 1
                break
            t = round(row[0], 9)
            if t not in ref[varnames[0]]:
                continue
            for iv, name in enumerate(varnames):
                for a, b in zip(row[1 + iv * nz:1 + (iv + 1) * nz], ref[name][t]):
                    if abs(a - b) > ATOL + RTOL * abs(b):
                        print("%s, t = %g, %s: %g in the log, %g from do_line_sampling"
                              % (log, t, name, a, b))
                        failed += 1
            nchecked += 1
        print("%s: %d rows of %d values (%d conserved components), %d compared with %s.*"
              % (log, len(rows), nvals, ncons, nchecked, prefix))
        if nchecked == 0:
            print("%s: no row has a matching do_line_sampling time" % log)
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
