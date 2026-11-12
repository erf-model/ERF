#!/usr/bin/env python3
"""Check a run of the example dust model (Tests/DustModelExample).

usage: check_dust_example.py run.log dust_example.csv FLUX AREA

Nothing leaves the domain of inputs_dust_example, so at every sum the dust in
the domain ERF prints (" RHO DUST   = M" after "TIME= t") must be the mass
the model emitted, FLUX * AREA * t, and the model's own count of it (the
CSV's emitted column, at the CSV row of time t) must agree.
"""

import re
import sys

REL_TOL = 1.0e-10


def sums(log):
    """(time, dust mass) for every sum ERF printed."""
    out, time = [], None
    with open(log) as f:
        for line in f:
            m = re.match(r"TIME=\s*([0-9.eE+-]+)", line)
            if m:
                time = float(m.group(1))
                continue
            m = re.match(r"\s*RHO DUST\s*=\s*([0-9.eE+-]+)", line)
            if m and time is not None:
                out.append((time, float(m.group(1))))
    return out


def emitted(csv):
    """time -> emitted mass, from the model's CSV."""
    rows = {}
    with open(csv) as f:
        next(f)
        for line in f:
            _, t, e = line.strip().split(",")
            rows[float(t)] = float(e)
    return rows


def main():
    log, csv, flux, area = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    s = sums(log)
    rows = emitted(csv)
    if len(s) < 2:
        print(f"FAIL: {len(s)} RHO DUST sums in {log}, need at least 2")
        return 1
    ok = True
    for t, m in s:
        expect = flux * area * t
        scale = max(abs(expect), flux * area)
        if abs(m - expect) > REL_TOL * scale:
            print(f"FAIL: t = {t}: RHO DUST = {m!r}, emitted {expect!r}")
            ok = False
        matched = [e for tt, e in rows.items() if abs(tt - t) <= 1.0e-9 * max(1.0, abs(t))]
        if matched and abs(matched[0] - expect) > REL_TOL * scale:
            print(f"FAIL: t = {t}: the model counts {matched[0]!r} emitted, expected {expect!r}")
            ok = False
    print(f"{len(s)} sums checked: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
