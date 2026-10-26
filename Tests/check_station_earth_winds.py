#!/usr/bin/env python3
"""Check the earth-relative wind columns of a station file against its grid winds.

Usage: check_station_earth_winds.py <station.dat> [--alpha DEG] [--rtol R]

The station must request x_velocity, y_velocity, u_east, v_north, wind_speed and
wind_dir at the same heights. With the grid x axis at --alpha degrees
counterclockwise from east (0 for a run without latitude/longitude arrays), every
row must satisfy

    u_east  = u cos(alpha) - v sin(alpha)
    v_north = u sin(alpha) + v cos(alpha)
    wind_speed = sqrt(u_east^2 + v_north^2)
    wind_dir   = 270 - atan2(v_north, u_east) in degrees, in [0, 360)

to --rtol relative to the wind speed (the file prints ten significant digits).
The header must report the same alpha. The check fails if the grid wind has no
y component anywhere, since then a swapped or unrotated v could pass.
Exit status 0 on pass, 1 on fail.
"""

import argparse
import math
import re
import sys

COLUMN = re.compile(r"^# column (\d+): (\S+)")
HEIGHT = re.compile(r"height=([-+0-9.eE]+)")
ALPHA = re.compile(r"grid x axis at ([-+0-9.eE]+) deg")
NAMES = ("x_velocity", "y_velocity", "u_east", "v_north", "wind_speed", "wind_dir")


def read_station(path):
    cols = {}      # (name, height) -> 0-based column
    alphas = []
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                m = COLUMN.match(line)
                if m:
                    h = HEIGHT.search(line)
                    cols[(m.group(2), float(h.group(1)) if h else None)] = int(m.group(1)) - 1
                    a = ALPHA.search(line)
                    if a:
                        alphas.append(float(a.group(1)))
                continue
            if line.strip():
                rows.append([float(x) for x in line.split()])
    return cols, alphas, rows


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("station_file")
    ap.add_argument("--alpha", type=float, default=0.0)
    ap.add_argument("--rtol", type=float, default=1.0e-8)
    a = ap.parse_args(argv)

    cols, alphas, rows = read_station(a.station_file)
    heights = sorted({h for (n, h) in cols if n == "u_east"})
    ok = True
    if not rows or not heights:
        print("FAIL: no rows or no u_east columns in %s" % a.station_file)
        return 1
    for n in NAMES:
        for h in heights:
            if (n, h) not in cols:
                print("FAIL: no %s column at height %g" % (n, h))
                return 1
    if not alphas or any(abs(x - a.alpha) > 1.0e-6 for x in alphas):
        print("FAIL: header alpha %s, expected %g" % (alphas, a.alpha))
        ok = False

    ca, sa = math.cos(math.radians(a.alpha)), math.sin(math.radians(a.alpha))
    worst = 0.0
    vmax = 0.0
    for r in rows:
        for h in heights:
            u, v = r[cols[("x_velocity", h)]], r[cols[("y_velocity", h)]]
            u_e, v_n = u * ca - v * sa, u * sa + v * ca
            speed = math.hypot(u_e, v_n)
            wdir = (270.0 - math.degrees(math.atan2(v_n, u_e))) % 360.0
            vmax = max(vmax, abs(v))
            scale = max(speed, 1.0e-30)
            got_dir = r[cols[("wind_dir", h)]]
            ddir = abs(got_dir - wdir)
            ddir = min(ddir, 360.0 - ddir)
            errs = [abs(r[cols[("u_east", h)]] - u_e) / scale,
                    abs(r[cols[("v_north", h)]] - v_n) / scale,
                    abs(r[cols[("wind_speed", h)]] - speed) / scale,
                    math.radians(ddir)]
            worst = max(worst, max(errs))
    print("rows %d, heights %s, alpha %g deg, max |y_velocity| %.4g, worst rel. error %.3e"
          % (len(rows), heights, a.alpha, vmax, worst))
    if vmax <= 0.0:
        print("FAIL: y_velocity is zero everywhere; the test cannot see v")
        ok = False
    if worst > a.rtol:
        print("FAIL: worst relative error %.3e > %.1e" % (worst, a.rtol))
        ok = False
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
