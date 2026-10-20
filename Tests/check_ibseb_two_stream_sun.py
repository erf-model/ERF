#!/usr/bin/env python3
"""Check that the building faces see the two-stream sun in a run.

    check_ibseb_two_stream_sun.py <run directory> --dt <level-0 step> --tau <sw_transmission>

Reads ibseb_buildings.csv (the faces' sun of every report) and
radiation_diag.csv (the two-stream sweep's domain-mean top-of-atmosphere
shortwave, S0 cos z) of a run of Tests/test_files/IBSEB_TwoStreamSun with
erf.ibseb.sun_mode = two_stream, and asserts, naming the defect each guards:
  1. the run reported more than two steps and the sweep ran at the matching
     times (nothing to compare is a failure, not a pass);
  2. the faces' top-of-atmosphere shortwave, dni / tau^(1/cos z) x cos z,
     equals the sweep's SW_TOA to 1e-5 at the start of every reported step (a
     row carries the end-of-step time and the sun the step used, from its
     start): a longitude passed in degrees, a missing start_datetime offset or
     the prescribed Spencer sun (1.5 degrees off on this date) all fail it;
  3. the sun is east of the meridian (azimuth in (0, 180) degrees): the run is
     at 08:20 local solar time, so a mirrored azimuth fails it.
Stdlib only.
"""
import argparse
import math
import sys


def read(path):
    with open(path) as f:
        hdr = f.readline().strip().split(",")
        return [dict(zip(hdr, line.strip().split(","))) for line in f if line.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--dt", type=float, required=True)
    ap.add_argument("--tau", type=float, required=True)
    a = ap.parse_args()
    ok = True

    def check(name, cond, detail):
        nonlocal ok
        print(f"  {name}: {'PASS' if cond else 'FAIL'} ({detail})")
        ok &= bool(cond)

    faces = [r for r in read(a.run + "/ibseb_buildings.csv") if r["level"] == "0" and r["building"] == "1"]
    sweep = {round(float(r["time"]), 6): float(r["SW_TOA"]) for r in read(a.run + "/radiation_diag.csv")
             if r["level"].strip() == "0"}
    pairs = []
    for r in faces:
        t = float(r["time_s"])
        cz = math.cos(math.radians(float(r["sun_zenith_deg"])))
        key = round(t - a.dt, 6)
        if t <= 0.0 or key not in sweep or cz <= 0.05:
            continue
        toa_faces = float(r["dni_Wm2"]) / a.tau ** (1.0 / cz) * cz
        pairs.append((t, toa_faces, sweep[key], float(r["sun_azimuth_deg"])))
    check("1. reports to compare", len(pairs) > 2, f"{len(pairs)} steps")
    worst = max((abs(p[1] / p[2] - 1.0) for p in pairs), default=float("inf"))
    check("2. the faces' S0 cos z is the sweep's", worst < 1e-5, f"largest relative difference {worst:.1e}")
    az = [p[3] for p in pairs]
    check("3. a morning sun in the east", bool(az) and all(0.0 < v < 180.0 for v in az),
          f"azimuth {min(az, default=float('nan')):.2f}-{max(az, default=float('nan')):.2f} deg")
    print("ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
