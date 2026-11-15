#!/usr/bin/env python3
"""Compare what the wrfinput reader reports for the land and soil fields with the file.

  check_log.py full LOG              every field read: components, min and max
  check_log.py partial LOG NAME...   the named fields skipped, the others read
  check_log.py off LOG               no field asked for
  check_log.py zero LOG              ZS and DZS all zero: taken as not given, the others read
  check_log.py bad LOG               the run stopped on ALBEDO
"""

import re
import sys

from make_wrfinput import DZS, NSOIL, ZS, surface_fields

RTOL = 2.0e-5   # the file holds 4-byte reals and the log prints 6 digits

FIELD = re.compile(r"wrfinput surface field (\w+) at level 0: (\d+) components?, "
                   r"min (\S+), max (\S+) \[")


def close(a, b):
    return abs(a - b) <= RTOL * max(abs(a), abs(b), 1.0e-3)


def reported(log):
    return {m.group(1): (int(m.group(2)), float(m.group(3)), float(m.group(4)))
            for m in FIELD.finditer(log)}


def check_fields(log, omitted, zero_layers=False):
    got = reported(log)
    errors = []
    for name, data in surface_fields().items():
        if name in omitted:
            if name in got:
                errors.append(f"{name}: reported although it is not in the file")
            if not re.search(rf"Skipping {name}\s*$", log, re.M):
                errors.append(f"{name}: no 'Skipping {name}' line")
            continue
        if name not in got:
            errors.append(f"{name}: not reported")
            continue
        ncomp, vmin, vmax = got[name]
        want = (NSOIL if data.ndim == 3 else 1, float(data.min()), float(data.max()))
        if ncomp != want[0] or not close(vmin, want[1]) or not close(vmax, want[2]):
            errors.append(f"{name}: reported {got[name]}, the file holds {want}")
    for label, name, values in (("depths ZS", "ZS", ZS), ("thicknesses DZS", "DZS", DZS)):
        match = re.search(rf"wrfinput soil layer {label} at level 0:(.*) \[m\]", log)
        if name in omitted:
            if match:
                errors.append(f"{name}: reported although it is not in the file")
            if not re.search(rf"Skipping {name}\s*$", log, re.M):
                errors.append(f"{name}: no 'Skipping {name}' line")
            continue
        if zero_layers:
            if match:
                errors.append(f"{name}: reported although it is zero in every layer")
            if f"wrfinput {name} at level 0 is zero in every soil layer: taken as not given" not in log:
                errors.append(f"{name}: no 'taken as not given' line")
            continue
        if not match:
            errors.append(f"{name}: not reported")
            continue
        got_layers = [float(v) for v in match.group(1).split()]
        if len(got_layers) != len(values) or not all(close(a, float(b)) for a, b in zip(got_layers, values)):
            errors.append(f"{name}: reported {got_layers}, the file holds {list(values)}")
    return errors


def main():
    mode, log = sys.argv[1], open(sys.argv[2]).read()
    if mode == "full":
        errors = check_fields(log, [])
    elif mode == "partial":
        errors = check_fields(log, sys.argv[3:])
    elif mode == "zero":
        errors = check_fields(log, [], zero_layers=True)
    elif mode == "off":
        errors = [f"{name}: asked for without the force-restore balance"
                  for name in list(surface_fields()) + ["ZS", "DZS"]
                  if f"Checking for {name} " in log]
        if "Checking for LANDMASK " not in log:
            errors.append("the reader's own field list was not found in the log")
    elif mode == "bad":
        errors = [] if "wrfinput ALBEDO holds a value outside its allowed range" in log else \
            ["the run did not stop on ALBEDO"]
        match = re.search(r"wrfinput ALBEDO at level 0: (\d+) value\(s\) outside 0 to 1 \[1\], (\d+) not finite; "
                          r"the finite values span (\S+) to (\S+)", log)
        if not match or (match.group(1), match.group(2)) != ("1", "0") or \
                not close(float(match.group(3)), 0.1) or not close(float(match.group(4)), 1.5):
            errors.append("the range report for ALBEDO is missing or wrong")
    else:
        raise SystemExit(f"unknown mode {mode}")
    for error in errors:
        print("FAIL:", error)
    if errors:
        raise SystemExit(1)
    print(f"PASS: {mode}")


if __name__ == "__main__":
    main()
