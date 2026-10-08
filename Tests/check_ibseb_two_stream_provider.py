#!/usr/bin/env python3
"""Check the building faces on the two-stream columns (erf.ibseb.radiation = two_stream).

    check_ibseb_two_stream_provider.py <work dir> --nz <cells in z> --cosz <cos zenith>
        --tau <SW optical depth per layer> --two-level-cosz <cos zenith>

Reads the per-face dumps (faces.stepNNNNNN.rank*.csv, faces/set[.lev1].rank*.csv) and the
per-building reports (ibseb_buildings.csv) of the runs Tests/RunIBSEBTwoStreamProvider.cmake
makes under <work dir>:
  - two_stream/: the cube of Tests/test_files/IBSEB_TwoStreamProvider on the two-stream
    columns under a transparent sky, two steps, then restarted from step 1 to step 2;
  - prescribed/: the prescribed provider set to the same sky (direct-normal irradiance S0, no
    diffuse light, no sky longwave, the column's ground albedo, emissivity and temperature);
  - absorbing/: the two_stream run with a shortwave optical depth of tau per layer and no
    scattering, so the beam at interface m is S0 cos z exp(-tau (nz - m) / cos z);
  - night/: IBSEB_TwoStreamProviderNight.i, the calendar sun below the horizon, traps on,
    an absorbing sky and a 70 m tower beside the 40 m cube;
  - two_level_clear/, two_level_absorbing/: Tests/test_files/IBSEB_RefinedLevels (a cube on
    level 1, a taller tower outside it) on the two-stream columns, transparent and absorbing.
A face takes its column at its own height: a roof the interface of its fluid cell's bottom
(m = k), a wall the mean of its cell's two (m = k and k + 1). It asserts, naming the
defect each guards:
  1. every run dumped the same faces, the cube's roof and four walls, with sunlit and
     shadowed faces and faces that see the ground (nothing to compare is a failure);
  2. on the transparent sky the two providers give every face the same direct and diffuse
     shortwave, absorbed shortwave and longwave at both steps (relative 1e-10): a beam not
     divided by cos z, the total down taken as the diffuse light, the ground's reflection or
     emission missing or taken at the wrong interface all fail it;
  3. under the absorbing sky every sunlit face's beam is the transparent one times
     exp(-tau (nz - m) / cos z) at its own sample, and the light the ground reflects onto the
     walls is times exp(-tau nz / cos z): one sample height for every face (the top of the
     buildings, the ground or the top of the atmosphere), or the ground's reflection taken
     above the ground, fails it;
  4. at night under the absorbing sky the roofs' sky longwave (LW_ext / f_sky, roofs see no
     ground) is positive and lower on the tower than on the cube, as less absorbing air lies
     above it: the longwave down read at one height for every face, or not at all, fails it;
  5. the initial report, before the first sweep, carries no radiation on the faces;
  6. on both levels of the two-level deck every sunlit face's beam follows check 3 and the
     ground's reflection too: a refined level sampling at a height of its own (the top of its
     own buildings, its own cell count) fails it. That a face reads its own column, not
     another one, the unit test IBSEBTwoStreamFaces checks on columns that all differ;
  7. at night the faces get no shortwave and the sky's longwave, and the run did not trap;
  8. a restart from step 1 neither appends a step-1 report with no radiation nor overwrites
     the step-1 face dump: the restart's initial report is not written with the two_stream
     provider, and the run goes on to report step 2.
Stdlib only.
"""
import argparse
import csv
import glob
import math
import sys


def faces(path_glob):
    """Faces of the dump files matching path_glob, keyed by (i, j, k, dir, side)."""
    rows = {}
    for path in sorted(glob.glob(path_glob)):
        with open(path) as f:
            hdr = f.readline().strip().split(",")
            for line in f:
                if not line.strip():
                    continue
                r = dict(zip(hdr, line.strip().split(",")))
                key = tuple(int(r[c]) for c in ("i", "j", "k", "dir", "side"))
                rows[key] = {c: float(v) for c, v in r.items()}
    return rows


def step_dump(run, step):
    return faces(f"{run}/faces.step{step:06d}.rank*.csv")


def close(a, b, rtol=1e-10, atol=1e-9):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def beam_factor(key, nz, tau, mu):
    """Transmission of the beam to a face's sample: a roof's interface k, a wall's mean of k, k+1."""
    k, d = key[2], key[3]
    t = lambda m: math.exp(-tau * (nz - m) / mu)
    return t(k) if d == 2 else 0.5 * (t(k) + t(k + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("--nz", type=int, required=True)
    ap.add_argument("--cosz", type=float, required=True)
    ap.add_argument("--tau", type=float, required=True)
    ap.add_argument("--two-level-cosz", type=float, required=True)
    a = ap.parse_args()
    ok = True

    def check(name, cond, detail):
        nonlocal ok
        print(f"  {name}: {'PASS' if cond else 'FAIL'} ({detail})")
        ok &= bool(cond)

    w = a.work
    ts, pre, ab = (step_dump(f"{w}/{r}", 1) for r in ("two_stream", "prescribed", "absorbing"))
    ts2, pre2 = step_dump(f"{w}/two_stream", 2), step_dump(f"{w}/prescribed", 2)
    sunlit = [k for k, f in ts.items() if f["SW_direct_in"] > 1.0]
    shaded = [k for k, f in ts.items() if f["shadow"] > 0.5]
    ground = [k for k, f in ts.items() if f["f_ground"] > 0.05]
    dirs = sorted({k[3] for k in ts})
    check("1. faces to compare", len(ts) > 0 and set(ts) == set(pre) == set(ab) == set(ts2) == set(pre2)
          and dirs == [0, 1, 2] and sunlit and shaded and ground,
          f"{len(ts)} faces, directions {dirs}, {len(sunlit)} sunlit, {len(shaded)} shadowed, "
          f"{len(ground)} seeing the ground")
    if not ok:
        # The checks below compare the runs face by face; without the same faces in each
        # there is nothing to compare.
        print("FAILED")
        return 1

    worst = None
    for x, y in ((ts, pre), (ts2, pre2)):
        for key, f in x.items():
            for c in ("SW_direct_in", "SW_diffuse_in", "SW_abs", "LW_in", "LW_net", "LW_ext"):
                if worst is None and not close(f[c], y[key][c]):
                    worst = (key, c, f[c], y[key][c])
    check("2. transparent sky: two_stream = prescribed", worst is None,
          "every face and term within 1e-10 at steps 1 and 2" if worst is None else f"first mismatch {worst}")

    refl = math.exp(-a.tau * a.nz / a.cosz)
    bad_dir = [k for k in sunlit
               if not close(ab[k]["SW_direct_in"], ts[k]["SW_direct_in"] * beam_factor(k, a.nz, a.tau, a.cosz), 1e-9)]
    bad_dif = [k for k in ts if not close(ab[k]["SW_diffuse_in"], ts[k]["SW_diffuse_in"] * refl, 1e-9)]
    heights = sorted({k[2] for k in sunlit})
    check("3. absorbing sky: the beam at each face's height, the reflection from the ground",
          len(heights) > 2 and not bad_dir and not bad_dif,
          f"{len(sunlit)} sunlit faces in {len(heights)} cells of height, reflection factor {refl:.6f}; "
          f"{len(bad_dir)} direct and {len(bad_dif)} diffuse mismatches"
          + (f", first direct {bad_dir[0]}: {ab[bad_dir[0]]['SW_direct_in']} vs "
             f"{ts[bad_dir[0]]['SW_direct_in'] * beam_factor(bad_dir[0], a.nz, a.tau, a.cosz)}" if bad_dir else ""))

    levels = {0: "faces/set", 1: "faces/set.lev1"}
    clear = {lev: faces(f"{w}/two_level_clear/{p}.rank*.csv") for lev, p in levels.items()}
    absd = {lev: faces(f"{w}/two_level_absorbing/{p}.rank*.csv") for lev, p in levels.items()}
    night = step_dump(f"{w}/night", 2)
    roof_lw = {}
    for f in night.values():
        if f["dir"] == 2 and f["f_sky"] > 0.5:
            roof_lw.setdefault(round(f["z_m"], 3), []).append(f["LW_ext"] / f["f_sky"])
    zs = sorted(roof_lw)
    means = [sum(roof_lw[z]) / len(roof_lw[z]) for z in zs]
    check("4. night: the roofs' sky longwave falls with height",
          len(zs) >= 2 and min(means) > 0.0 and all(m1 > m2 for m1, m2 in zip(means, means[1:])),
          ", ".join(f"{z:g} m: {m:.2f} W/m2" for z, m in zip(zs, means)))

    init = step_dump(f"{w}/two_stream", 0)
    lit0 = [k for k, f in init.items() if f["SW_abs"] != 0.0 or f["LW_in"] != 0.0]
    check("5. no radiation before the first sweep", init and not lit0,
          f"{len(init)} faces at step 0, {len(lit0)} with radiation")

    mu = a.two_level_cosz
    detail, good = [], bool(clear[0]) and bool(clear[1])
    for lev in levels:
        c, b = clear[lev], absd[lev]
        if set(c) != set(b) or not c:
            good = False
            detail.append(f"level {lev}: {len(c)} and {len(b)} faces in the two runs")
            continue
        lit = [k for k, f in c.items() if f["SW_direct_in"] > 1.0]
        bad = [k for k in lit
               if not close(b[k]["SW_direct_in"], c[k]["SW_direct_in"] * beam_factor(k, a.nz, a.tau, mu), 1e-9)]
        bad += [k for k in c
                if not close(b[k]["SW_diffuse_in"], c[k]["SW_diffuse_in"] * math.exp(-a.tau * a.nz / mu), 1e-9)]
        good &= bool(lit) and not bad
        detail.append(f"level {lev}: {len(c)} faces, {len(lit)} sunlit, {len(bad)} mismatches")
    check("6. two levels: each level's columns at each face's height", good, "; ".join(detail))

    lit_n = [k for k, f in night.items() if f["SW_abs"] != 0.0]
    cold_n = [k for k, f in night.items() if f["dir"] == 2 and not f["LW_ext"] > 0.0]
    check("7. night: no shortwave, the sky's longwave", night and not lit_n and not cold_n,
          f"{len(night)} faces at step 2, {len(lit_n)} with shortwave, {len(cold_n)} roofs without sky longwave")

    with open(f"{w}/two_stream/ibseb_buildings.csv") as f:
        rows = list(csv.DictReader(f))
    s1 = [r for r in rows if r["step"] == "1"]
    s2 = [r for r in rows if r["step"] == "2"]
    dump1 = step_dump(f"{w}/two_stream", 1)
    check("8. a restart keeps the step-1 report and dump",
          len(s1) == 1 and float(s1[0]["SW_abs_mean_Wm2"]) > 0.0 and len(s2) == 2
          and any(f["SW_abs"] > 0.0 for f in dump1.values()),
          f"{len(s1)} step-1 row(s), SW_abs_mean {[r['SW_abs_mean_Wm2'] for r in s1]}; "
          f"{len(s2)} step-2 rows (before and after the restart); step-1 dump SW_abs max "
          f"{max((f['SW_abs'] for f in dump1.values()), default=float('nan')):.1f}")

    print("ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
