#!/usr/bin/env python3
"""Compare what the wrfinput reader reports for the land and soil fields with the file.

  check_log.py full LOG              every field read: components, min and max
  check_log.py partial LOG NAME...   the named fields skipped, the others read
  check_log.py off LOG               no field asked for
  check_log.py zero LOG              ZS and DZS all zero: taken as not given, the others read
  check_log.py surfacezs LOG         ZS starts at zero: held
  check_log.py drop LOG NAME=N...    the named fields reported with N values outside their
                                     range and not held; the others read; the run completed
  check_log.py restart LOG           a restart: every field restored from the checkpoint
  check_log.py regrid LOG            level 1 made by a regrid, then moved and widened by the
                                     next: every field of level 0 on its grids both times
  check_log.py tagged LOG            level 1 made by tagging at the start: every field of
                                     level 0 on its grids
"""

import re
import sys

from make_wrfinput import DZS, NSOIL, ZS, ZS_SURFACE, surface_fields

RTOL = 2.0e-5   # the file holds 4-byte reals and the log prints 6 digits

FIELD = re.compile(r"wrfinput surface field (\w+) at level 0: (\d+) components?, "
                   r"min (\S+), max (\S+) \[")
RESTORED = re.compile(r"wrfinput surface field (\w+) at level 0 restored from the checkpoint: "
                      r"(\d+) components?, min (\S+), max (\S+) \[")


ON_NEW_GRIDS = re.compile(r"wrfinput surface field (\w+) at level 1 on the new grids: "
                          r"(\d+) components?, (\d+) of (\d+) columns from level 0, "
                          r"min (\S+), max (\S+) \[")

# Level 1 (refined by 2 in x and y) in the regrid leg: the level 0 columns under it as
# (j range, i range), the fine columns taken from level 0 and the fine columns in all.
# First the level is new, on level 0 columns 2..5 x 2..5: all 64 from level 0. Then it
# moves to 4..9 x 2..5: its old grids give 32 of the 96, level 0 the other 64.
REGRID_STAGES = (((2, 6), (2, 6), 64, 64), ((2, 6), (4, 10), 64, 96))


def close(a, b):
    return abs(a - b) <= RTOL * max(abs(a), abs(b), 1.0e-3)


def reported(log, pattern=FIELD):
    return {m.group(1): (int(m.group(2)), float(m.group(3)), float(m.group(4)))
            for m in pattern.finditer(log)}


def check_fields(log, omitted, zero_layers=False, dropped=None, depths=ZS, restored=False):
    got = reported(log, RESTORED if restored else FIELD)
    dropped = dropped or {}
    tail = " restored from the checkpoint" if restored else ""
    errors = []
    for name, data in surface_fields().items():
        if name in dropped:
            if name in got:
                errors.append(f"{name}: reported as held although it holds a value outside its range")
            if not re.search(rf"Warning: wrfinput {name} at level 0 is not held: {dropped[name]} "
                             rf"value\(s\) outside ", log):
                errors.append(f"{name}: no 'is not held' line with {dropped[name]} value(s) outside")
            continue
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
    for label, name, values in (("depths ZS", "ZS", depths), ("thicknesses DZS", "DZS", DZS)):
        match = re.search(rf"wrfinput soil layer {label} at level 0{tail}:(.*) \[m\]", log)
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


def check_regrid(log, stages=REGRID_STAGES, last_step=3):
    errors = []
    got = {}
    for m in ON_NEW_GRIDS.finditer(log):
        got.setdefault(m.group(1), []).append(
            (int(m.group(2)), int(m.group(3)), int(m.group(4)), float(m.group(5)), float(m.group(6))))
    for name, data in surface_fields().items():
        lines = got.get(name, [])
        if len(lines) != len(stages):
            errors.append(f"{name}: {len(lines)} report(s) at level 1, expected {len(stages)}")
            continue
        for stage, line in zip(stages, lines):
            (j0, j1), (i0, i1), from_coarse, columns = stage
            under = data[..., j0:j1, i0:i1]
            want = (NSOIL if data.ndim == 3 else 1, from_coarse, columns)
            if line[:3] != want or not close(line[3], float(under.min())) or \
                    not close(line[4], float(under.max())):
                errors.append(f"{name}: reported {line}, expected {want} and the range "
                              f"{float(under.min())} to {float(under.max())}")
    if re.search(r"Warning: wrfinput \w+ at level \d+ is dropped", log):
        errors.append("a field was dropped on a regrid")
    if last_step > 0:
        if f"Coarse STEP {last_step} ends." not in log:
            errors.append(f"the run did not reach the end of step {last_step}")
        if "Done with advance_dycore at level 1" not in log:
            errors.append("level 1 took no step")
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
    elif mode == "surfacezs":
        errors = check_fields(log, [], depths=ZS_SURFACE)
    elif mode == "drop":
        dropped = dict(item.split("=") for item in sys.argv[3:])
        errors = check_fields(log, [], dropped=dropped)
        if "ALBEDO" in dropped:
            match = re.search(r"wrfinput ALBEDO at level 0 is not held: (\d+) value\(s\) outside 0 to 1 \[1\], "
                              r"(\d+) not finite; the finite values span (\S+) to (\S+)", log)
            if not match or (match.group(1), match.group(2)) != ("1", "0") or \
                    not close(float(match.group(3)), 0.1) or not close(float(match.group(4)), 1.5):
                errors.append("the range report for ALBEDO is missing or wrong")
    elif mode == "restart":
        errors = check_fields(log, [], restored=True)
        if "Checking for ALBEDO " in log:
            errors.append("the restart read the fields from the wrfinput file again")
    elif mode == "regrid":
        errors = check_fields(log, []) + check_regrid(log)
    elif mode == "tagged":
        errors = check_fields(log, []) + check_regrid(log, REGRID_STAGES[:1], last_step=0)
    else:
        raise SystemExit(f"unknown mode {mode}")
    for error in errors:
        print("FAIL:", error)
    if errors:
        raise SystemExit(1)
    print(f"PASS: {mode}")


if __name__ == "__main__":
    main()
