#!/usr/bin/env python3
"""Compare two face dumps of the FaceStorage deck face by face.

    python3 compare_faces.py faces_straight_prognostic faces_restart
    python3 compare_faces.py faces_restart faces_restart_np1 --rtol 1e-8 --resid-atol 1e-6

Reads every ``<prefix>.rank<N>.csv`` of each prefix (one file per rank, so a
run on any number of ranks gives the same set of rows), keys the rows by
(i, j, k, dir, side) and requires the two sets to hold the same faces with
every column equal to the 12 significant digits the dump is written with.
The checkpoint stores the state in double precision and the step after the
restart sees the same atmosphere, so a difference means a record was lost,
misplaced or rounded. Exit code 0 on a
match, 1 otherwise, printing the first difference of each column.

With ``--rtol R`` a column passes when |x - y| <= R max(|x|, |y|), and
``--resid-atol A`` gives the balance residual ``resid`` an absolute band
instead: it is |balance| after the Newton solve, a cancellation of terms of
hundreds of W/m2 down to about 1e-8, so it carries no significant digits to
compare. The two are for runs on different rank counts, whose atmospheres
differ in the last digit (the ERF solve itself, not the face state); a lost or
misplaced record moves the skin temperature by kelvins and still fails.
"""
import argparse, csv, glob, sys

def read_dump(prefix):
    files = sorted(glob.glob(f"{prefix}.rank*.csv"))
    if not files:
        sys.exit(f"no files {prefix}.rank*.csv")
    rows, header = {}, None
    for fn in files:
        with open(fn) as f:
            r = csv.reader(f)
            h = next(r)
            if header is None:
                header = h
            elif h != header:
                sys.exit(f"{fn}: header differs from {files[0]}")
            for row in r:
                key = tuple(int(v) for v in row[:5])
                if key in rows:
                    sys.exit(f"{fn}: face {key} appears twice")
                rows[key] = row
    return header, rows

def equal(name, x, y, rtol, resid_atol):
    if x == y:
        return True
    fx, fy = float(x), float(y)
    if name == "resid" and resid_atol > 0:
        return abs(fx - fy) <= resid_atol
    return abs(fx - fy) <= rtol * max(abs(fx), abs(fy))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("a"); ap.add_argument("b")
    ap.add_argument("--rtol", type=float, default=0.0, help="relative tolerance per column (default exact)")
    ap.add_argument("--resid-atol", type=float, default=0.0, help="absolute tolerance for the resid column (default exact)")
    args = ap.parse_args()
    a, b = args.a, args.b
    ha, ra = read_dump(a)
    hb, rb = read_dump(b)
    ok = True
    if ha != hb:
        print("headers differ"); ok = False
    if set(ra) != set(rb):
        only_a, only_b = sorted(set(ra) - set(rb)), sorted(set(rb) - set(ra))
        print(f"face sets differ: {len(only_a)} only in {a} (first {only_a[:3]}), {len(only_b)} only in {b} (first {only_b[:3]})")
        ok = False
    bad = {}
    for key in sorted(set(ra) & set(rb)):
        for name, x, y in zip(ha, ra[key], rb[key]):
            if name not in bad and not equal(name, x, y, args.rtol, args.resid_atol):
                bad[name] = (key, x, y)
    for name, (key, x, y) in bad.items():
        print(f" {name}: face {key} {x} vs {y}"); ok = False
    how = "all columns equal" if args.rtol == 0 and args.resid_atol == 0 else f"all columns equal to rtol {args.rtol:g} (resid within {args.resid_atol:g})"
    print(f"{len(ra)} faces in {a}, {len(rb)} in {b}: " + (how if ok else "DIFFER"))
    sys.exit(0 if ok else 1)
