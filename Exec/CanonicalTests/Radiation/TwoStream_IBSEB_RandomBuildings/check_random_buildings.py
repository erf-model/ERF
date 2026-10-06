#!/usr/bin/env python3
"""Check the random-building morning under the two-stream radiation, on two levels.

    python3 check_random_buildings.py [run directory]       # default: the current directory

Reads buildings.csv (make_buildings.py), ibseb_buildings.csv (one row per
building and level every erf.ibseb.csv_int steps), the face dumps
faces/set.step* (level 0) and faces/set.lev1.step* (level 1),
radiation_diag.csv (the two-stream diagnostics) and the 2D plotfiles plt2d*.
Each check names the defect it guards; the exit code is the number of failed
checks.

  1. both levels hold every building of buildings.csv, numbered alike: a
     building wholly inside the refined level keeps its own faces there, and
     the ids follow the same scan order;
  2. the balance closes on every building and level all morning (residual
     below 1e-3 W/m2);
  3. the faces see the two-stream sun, east of the meridian all morning
     (azimuth in (0, 180) degrees, so a mirrored sun fails): the top-of-atmosphere irradiance the
     faces' direct beam implies (dni / tau^(1/cos z)) equals the two-stream
     sweep's SW_TOA / cos z, to 1e-5. A CSV row is stamped with the end of its
     step and carries the sun the step used, from its start, so the sweep's
     row one step (erf.fixed_dt) earlier is the one compared: the same sun
     one step late already misses by up to 7e-5, and the faces' own prescribed
     sun (sun_mode = solar, with the equation of time) by 2.4 % at 08:20;
  4. the sun climbs: the zenith falls through the run and the buildings'
     mean shadow fraction ends below where it starts, on both levels;
  5. the buildings warm: every building's mean skin temperature ends above
     the 300 K it started at, on both levels;
  6. the two levels agree: each building's final mean skin temperature
     within 3 K between the levels (the refined level resolves the rims and
     the air next to the walls better). The face areas are reported, not
     checked: on level 0 a building also fills the cell its roof ramps down
     over, which adds most to the smallest buildings;
  7. the ground heats outside the buildings: over the refined area, away from
     the building footprints and the cell around them, the force-restore skin
     of level 1 ends above 300 K, and so does level 0's outside the refined
     area. The two means are reported, not compared: under the refined area
     level 0 holds level 1's average (the average down), and the open ground
     outside it differs by the buildings' influence. The two-stream ground
     balance still covers the footprints; they are left out here and reported
     separately with their sensible heat flux.
"""
import glob
import os
import re
import sys

import numpy as np


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def read_csv(path):
    a = np.genfromtxt(path, delimiter=",", names=True)
    return {n: a[n] for n in a.dtype.names}


def dump_steps(prefix):
    """Steps of the tagged face dumps of one level prefix (faces/set or faces/set.lev1)."""
    steps = set()
    for fn in glob.glob(prefix + ".step*.rank*.csv"):
        m = re.search(re.escape(prefix) + r"\.step(\d+)\.rank\d+\.csv$", fn)
        if m:
            steps.add(int(m.group(1)))
    return sorted(steps)


def read_dump(prefix, step):
    rows, hdr = [], None
    for fn in sorted(glob.glob(f"{prefix}.step{step:06d}.rank*.csv")):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            rows += [[float(v) for v in line.split(",")] for line in f if line.strip()]
    a = np.array(rows)
    return {h: a[:, n] for n, h in enumerate(hdr)}


def read_plotfile_2d(path):
    """Fields of an AMReX plotfile written by ERF's plt2d, per level, as
    (names, [(lo_i, lo_j, array[ncomp, ni, nj]) per box], dx) per level."""
    with open(os.path.join(path, "Header")) as f:
        lines = f.read().split("\n")
    nvar = int(lines[1])
    names = lines[2:2 + nvar]
    p = 2 + nvar
    p += 1                          # space dimension
    p += 1                          # time
    finest = int(lines[p]); p += 1
    p += 3                          # prob_lo, prob_hi, ref ratios
    p += 1                          # domains
    p += 1                          # level steps
    dxs = []
    for lev in range(finest + 1):
        dxs.append([float(v) for v in lines[p].split()]); p += 1
    out = []
    for lev in range(finest + 1):
        boxes = []
        with open(os.path.join(path, f"Level_{lev}", "Cell_H")) as f:
            ch = f.read()
        for name, off in re.findall(r"FabOnDisk: (\S+) (\d+)", ch):
            with open(os.path.join(path, f"Level_{lev}", name), "rb") as f:
                f.seek(int(off))
                head = b""
                while not head.endswith(b"\n"):
                    head += f.read(1)
                m = re.search(rb"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\) \(\d+,\d+,\d+\)\) (\d+)", head)
                lo = [int(m.group(k)) for k in (1, 2, 3)]
                hi = [int(m.group(k)) for k in (4, 5, 6)]
                nc = int(m.group(7))
                ni, nj, nk = hi[0] - lo[0] + 1, hi[1] - lo[1] + 1, hi[2] - lo[2] + 1
                order = "<" if b"(8 7 6 5 4 3 2 1)" in head else ">"
                data = np.frombuffer(f.read(8 * ni * nj * nk * nc), dtype=order + "f8")
                arr = data.reshape(nc, nk, nj, ni)[:, 0, :, :].transpose(0, 2, 1)   # comp, i, j
                boxes.append((lo[0], lo[1], arr))
        out.append(boxes)
    return names, out, dxs


def level_field(boxes, comp, shape):
    """One component of a level's boxes on the level's index space, NaN where it has no cells."""
    f = np.full(shape, np.nan)
    for i0, j0, arr in boxes:
        f[i0:i0 + arr.shape[1], j0:j0 + arr.shape[2]] = arr[comp]
    return f


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def main():
    run = sys.argv[1] if len(sys.argv) > 1 else "."
    os.chdir(run)
    nfail = 0

    def check(name, ok, detail):
        nonlocal nfail
        print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})")
        nfail += 0 if ok else 1

    bld = read_csv("buildings.csv")
    nb = len(bld["building"])
    c = read_csv("ibseb_buildings.csv")
    lev = c["level"].astype(int)
    tau, dt = None, None
    with open("inputs") as f:
        for line in f:
            m = re.match(r"\s*erf\.ibseb\.sw_transmission\s*=\s*(\S+)", line)
            if m:
                tau = float(m.group(1))
            m = re.match(r"\s*erf\.fixed_dt\s*=\s*(\S+)", line)
            if m:
                dt = float(m.group(1))
    tau = 0.7 if tau is None else tau          # IBSEBParams default
    if dt is None:
        raise SystemExit("inputs: erf.fixed_dt is not set; check 3 needs the level-0 step")

    # 1. buildings on both levels, numbered alike
    s0, s1 = dump_steps("faces/set"), dump_steps("faces/set.lev1")
    d0, d1 = read_dump("faces/set", s0[0]), read_dump("faces/set.lev1", s1[0])
    xc = 0.5 * (bld["x_lo_m"] + bld["x_hi_m"]); yc = 0.5 * (bld["y_lo_m"] + bld["y_hi_m"])
    worst = 0.0
    for d in (d0, d1):
        for b in range(1, nb + 1):
            roof = (d["bid"] == b) & (d["dir"] == 2)
            if roof.any():
                worst = max(worst, abs(d["x_m"][roof].mean() - xc[b - 1]), abs(d["y_m"][roof].mean() - yc[b - 1]))
            else:
                worst = np.inf
    ok = int(d0["bid"].max()) == nb and int(d1["bid"].max()) == nb and worst <= 10.0
    check("1. every building on both levels, numbered alike", ok,
          f"{int(d0['bid'].max())} and {int(d1['bid'].max())} buildings of {nb}; roof centres within {worst:.1f} m of buildings.csv")

    # 2. residual
    check("2. the balance closes", c["resid_max_Wm2"].max() < 1e-3,
          f"largest residual {c['resid_max_Wm2'].max():.1e} W/m2 over {len(lev)} rows")

    # 3. the two-stream sun
    r = read_csv("radiation_diag.csv")
    rt = r["time"][r["level"] == 0]; rtoa = r["SW_TOA"][r["level"] == 0]
    first = (lev == 0) & (c["building"] == 1)
    ct, cz = c["time_s"][first], np.cos(np.radians(c["sun_zenith_deg"][first]))
    dni = c["dni_Wm2"][first]
    rel = []
    for t, z, q in zip(ct, cz, dni):
        # A CSV row carries the time at the end of its step and the sun the
        # step used, placed at the step's start: compare there.
        k = np.nonzero(np.abs(rt - (t - dt)) < 1e-6)[0]
        if len(k) == 0 or z <= 0.05:
            continue
        s0_faces = q / tau ** (1.0 / z)
        s0_cols = rtoa[k[0]] / z
        rel.append(abs(s0_faces / s0_cols - 1.0))
    rel = np.array(rel)
    az = c["sun_azimuth_deg"][first]
    check("3. the faces see the two-stream sun", len(rel) > 2 and rel.max() < 1e-5 and np.all((az > 0.0) & (az < 180.0)),
          f"{len(rel)} common times, largest relative difference {rel.max() if len(rel) else float('nan'):.1e}; "
          f"azimuth {az.min():.1f}-{az.max():.1f} deg")

    # 4. the sun climbs, the shadows shorten
    zen = c["sun_zenith_deg"][first]
    ok, det = zen[-1] < zen[0], [f"zenith {zen[0]:.1f} -> {zen[-1]:.1f} deg"]
    for L in (0, 1):
        m = lev == L
        t = c["time_s"][m]
        a = c["area_m2"][m]; sh = c["shadow_frac"][m]
        t0, t1 = t.min(), t.max()
        f0 = np.sum((sh * a)[t == t0]) / np.sum(a[t == t0]); f1 = np.sum((sh * a)[t == t1]) / np.sum(a[t == t1])
        ok &= f1 < f0
        det.append(f"level {L} shadow {f0:.3f} -> {f1:.3f}")
    check("4. the sun climbs and the shadows shorten", ok, ", ".join(det))

    # 5. the buildings warm
    ok, det = True, []
    for L in (0, 1):
        m = (lev == L) & (c["time_s"] == c["time_s"][lev == L].max())
        ok &= bool(np.all(c["T_skin_mean_K"][m] > 300.0))
        det.append(f"level {L} final means {c['T_skin_mean_K'][m].min():.1f}-{c['T_skin_mean_K'][m].max():.1f} K")
    check("5. every building warms", ok, ", ".join(det))

    # 6. the levels agree
    tend = min(c["time_s"][lev == 0].max(), c["time_s"][lev == 1].max())
    darea, dT = [], []
    for b in range(1, nb + 1):
        m0 = (lev == 0) & (c["building"] == b) & (c["time_s"] == tend)
        m1 = (lev == 1) & (c["building"] == b) & (c["time_s"] == tend)
        if not m0.any() or not m1.any():
            darea.append(np.inf); dT.append(np.inf)    # a building missing on a level: check 1's failure, reported here too
            continue
        darea.append(abs(c["area_m2"][m1][0] / c["area_m2"][m0][0] - 1.0))
        dT.append(c["T_skin_mean_K"][m1][0] - c["T_skin_mean_K"][m0][0])
    darea, dT = np.array(darea), np.array(dT)
    check("6. the two levels agree on every building", np.abs(dT).max() <= 3.0,
          f"final skin level 1 - level 0 from {dT.min():+.2f} to {dT.max():+.2f} K; face areas differ by up to {100 * darea.max():.0f} %")

    # 7. the ground away from the buildings
    plts = sorted(glob.glob("plt2d*"), key=lambda s: int(re.sub(r"\D", "", s)))
    names, levels, dxs = read_plotfile_2d(plts[-1])
    k = names.index("seb_t_sfc")
    nx0 = int(round(640.0 / dxs[0][0])); ny0 = int(round(640.0 / dxs[0][1]))
    f0 = level_field(levels[0], k, (nx0, ny0))
    f1 = level_field(levels[1], k, (2 * nx0, 2 * ny0))
    # Keep-out: every footprint grown by one coarse cell (its ramp) and one more.
    def keep_out(nx, ny, dx, dy):
        x = (np.arange(nx) + 0.5) * dx; y = (np.arange(ny) + 0.5) * dy
        m = np.zeros((nx, ny), dtype=bool)
        for b in range(nb):
            g = 2.0 * 20.0
            m |= np.outer((x > bld["x_lo_m"][b] - g) & (x < bld["x_hi_m"][b] + g),
                          (y > bld["y_lo_m"][b] - g) & (y < bld["y_hi_m"][b] + g))
        return m
    ko1 = keep_out(2 * nx0, 2 * ny0, dxs[1][0], dxs[1][1])
    ko0 = keep_out(nx0, ny0, dxs[0][0], dxs[0][1])
    on1 = ~np.isnan(f1)
    # The coarse columns under the refined level.
    on0 = np.zeros_like(ko0)
    for i0, j0, arr in levels[1]:
        on0[i0 // 2:(i0 + arr.shape[1]) // 2, j0 // 2:(j0 + arr.shape[2]) // 2] = True
    g1, g_out = f1[on1 & ~ko1].mean(), f0[~on0].mean()
    under = f1[on1 & ko1].mean()
    hf = names.index("seb_hfx") if "seb_hfx" in names else None
    flux = ""
    if hf is not None:
        h1 = level_field(levels[1], hf, (2 * nx0, 2 * ny0))
        flux = f"; sensible heat flux {np.nanmean(h1[on1 & ~ko1]):.0f} W/m2 there, {np.nanmean(h1[on1 & ko1]):.0f} W/m2 under and around the buildings"
    check("7. the open ground heats",
          g1 > 300.0 and g_out > 300.0,
          f"{plts[-1]}: open ground {g1:.2f} K in the refined area (level 1) and {g_out:.2f} K outside it (level 0); "
          f"the footprints and their surroundings, left out, {under:.2f} K on level 1{flux}")

    # Not a check: where each building's absorbed shortwave comes from at the
    # last dump, on both levels, the numbers the README's reading rests on.
    print("  breakdown at the last face dump (level 0 / level 1): roof share of the face area, "
          "absorbed SW on roofs and on walls [W/m2], sunlit share of the wall area")
    last = {0: read_dump("faces/set", s0[-1]), 1: read_dump("faces/set.lev1", s1[-1])}
    for b in range(1, nb + 1):
        cols = []
        for L in (0, 1):
            d = last[L]
            m = d["bid"] == b
            roof, wall = m & (d["dir"] == 2), m & (d["dir"] != 2)
            a = d["area_m2"]
            share = a[roof].sum() / a[m].sum()
            sw_r = (d["SW_abs"][roof] * a[roof]).sum() / a[roof].sum()
            sw_w = (d["SW_abs"][wall] * a[wall]).sum() / a[wall].sum()
            lit_w = (a[wall] * (d["SW_direct_in"][wall] > 0)).sum() / a[wall].sum()
            cols.append((share, sw_r, sw_w, lit_w))
        print(f"    building {b}: roof share {cols[0][0]:.2f} / {cols[1][0]:.2f}, "
              f"roof SW {cols[0][1]:.0f} / {cols[1][1]:.0f}, wall SW {cols[0][2]:.0f} / {cols[1][2]:.0f}, "
              f"sunlit walls {cols[0][3]:.2f} / {cols[1][3]:.2f}")

    print("ALL PASS" if nfail == 0 else f"{nfail} FAILED")
    return nfail


if __name__ == "__main__":
    sys.exit(main())
