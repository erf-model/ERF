#!/usr/bin/env python3
"""Check the building-set morning.

    python3 check_buildingset.py ibseb_set.csv faces/set run_set.log

From the per-building CSV (one row a minute), the face dumps (one every
five minutes) and the run log:
  1. four buildings with the expected face counts and materials
     (1 concrete slab, 2 timber block, 3 brick cube, 4 timber block);
  2. the balance residual stays below 1e-3 W/m2 on every building;
  3. mutual shadowing: at sunrise the slab is the most shaded building
     (the cube east of it throws its shadow onto the slab's east wall) and
     clears by late morning; the 20 m blocks are free of shadow by then;
  4. the building view fraction: the facing walls of the slab and the
     cube see more building than the far block's walls;
  5. materials show in the response: at the end the timber blocks' mean
     skin temperature exceeds the concrete slab's on the sunlit side
     (light cladding warms faster);
  6. the wall function beyond neutral is active: w* is positive on sunlit
     faces and most roofs are unstable by the end;
  7. the run log's timing line is reported (faces per rank, seconds per
     step of the balance), for estimating city-scale cost.
"""
import sys, glob, re
import numpy as np

def load(prefix):
    rows = []
    for fn in sorted(glob.glob(prefix + ".rank*.csv")):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            for line in f: rows.append([float(v) for v in line.strip().split(",")])
    a = np.array(rows); d = {h: a[:, n] for n, h in enumerate(hdr)}
    key = np.lexsort((d["side"], d["dir"], d["k"], d["j"], d["i"]))
    return {h: v[key] for h, v in d.items()}

def load_steps(prefix):
    files = sorted(set(re.sub(r"\.rank\d+\.csv$", "", fn) for fn in glob.glob(prefix + ".step*.rank*.csv")))
    return np.array([int(re.search(r"\.step(\d+)$", f).group(1)) for f in files]), files

def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})"); return ok

def main():
    csv, prefix, log = sys.argv[1:4]
    c = np.genfromtxt(csv, delimiter=",", names=True)
    steps, files = load_steps(prefix); th = 5.0 + (steps + 1) * 0.5 / 3600.0   # solar time [h], run starts at 05:00
    d0 = load(files[0]); dl = load(files[-1])
    ok = True
    # 1. buildings and materials
    nb = int(c["building"].max())
    counts = {b: int((d0["bid"] == b).sum()) for b in range(1, nb + 1)}
    mats = {b: int(np.unique(d0["mat"][d0["bid"] == b])[0]) for b in range(1, nb + 1)}
    ok &= report("four buildings with their materials", nb == 4 and mats == {1: 1, 2: 3, 3: 2, 4: 3},
                 ", ".join(f"building {b}: {counts[b]} faces, material {mats[b]}" for b in counts))
    # 2. residual
    ok &= report("balance residual on every building all morning", c["resid_max_Wm2"].max() < 1e-3, f"max {c['resid_max_Wm2'].max():.1e} W/m2 over {len(c)} rows")
    # 3. shadowing from the CSV shadow fraction per building. The sun rises in
    # the east-north-east, so at sunrise the 40 m cube throws its shadow west
    # onto the slab's east wall (the slab shades the cube in the afternoon,
    # outside this run). Every building also shades its own stepped rim at
    # low sun; the 20 m blocks lose that by mid-morning, the cube's 20 m step
    # keeps a little of it at 11:00.
    def col(b, name): return c[name][c["building"] == b]
    t = col(1, "time_s") / 3600.0 + 5.0
    sh = {b: col(b, "shadow_frac") for b in (1, 2, 3, 4)}
    early = (t > 5.4) & (t < 6.2); late = t > 10.4
    ok &= report("at sunrise the slab is the most shaded building (the cube's shadow on its east wall) and clears by late morning",
                 sh[1][early].mean() > max(sh[b][early].mean() for b in (2, 3, 4)) and sh[1][early].mean() - sh[1][late].mean() > 0.08,
                 f"slab {sh[1][early].mean():.3f} at 05:30-06:10 vs {sh[1][late].mean():.3f} after 10:30; others early " + ", ".join(f"{sh[b][early].mean():.3f}" for b in (2, 3, 4)))
    ok &= report("the 20 m blocks are free of shadow by late morning", sh[2][late].max() == 0.0 and sh[4][late].max() == 0.0,
                 f"north block {sh[2][late].max():.3f}, far block {sh[4][late].max():.3f} after 10:30")
    # 4. building view fraction
    slab_east = (d0["bid"] == 1) & (d0["dir"] == 0) & (d0["side"] == -1)
    cube_west = (d0["bid"] == 3) & (d0["dir"] == 0) & (d0["side"] == 1)
    far = (d0["bid"] == 4) & (d0["dir"] < 2)
    # Every wall sees some building (its own stepped rim, and on this periodic
    # domain the others across the boundary); the facing walls see clearly more.
    ok &= report("the facing walls of the slab and the cube see more building than the far block's walls",
                 d0["f_bldg"][slab_east].mean() > d0["f_bldg"][far].mean() + 0.05 and d0["f_bldg"][cube_west].mean() > d0["f_bldg"][far].mean() + 0.05,
                 f"slab east {d0['f_bldg'][slab_east].mean():.2f}, cube west {d0['f_bldg'][cube_west].mean():.2f}, far block walls {d0['f_bldg'][far].mean():.3f}")
    # 5. materials in the response: mean skin at the end, roofs only (all see the sun)
    roof = lambda b: (dl["bid"] == b) & (dl["dir"] == 2)
    ok &= report("timber roofs end warmer than the concrete roof (light cladding warms faster)",
                 dl["T_skin"][roof(2)].mean() > dl["T_skin"][roof(1)].mean() and dl["T_skin"][roof(4)].mean() > dl["T_skin"][roof(1)].mean(),
                 f"roof means at {th[-1]:.1f} h: slab {dl['T_skin'][roof(1)].mean():.1f}, north block {dl['T_skin'][roof(2)].mean():.1f}, cube {dl['T_skin'][roof(3)].mean():.1f}, far block {dl['T_skin'][roof(4)].mean():.1f} K")
    # 6. the wall function beyond neutral
    lit = dl["SW_direct_in"] > 100.0
    roofs = dl["dir"] == 2; unstable = (dl["Olen"][roofs] < 0).mean()
    ok &= report("w* positive on the sunlit faces and most roofs unstable by the end (the shaded rim roofs may be stable)",
                 (dl["w_star"][lit] > 0).all() and unstable > 0.8,
                 f"w* {dl['w_star'][lit].min():.2f}-{dl['w_star'][lit].max():.2f} m/s on {lit.sum()} sunlit faces, {100*unstable:.0f} % of {roofs.sum()} roofs unstable")
    # 7. timing
    tl = re.findall(r"\[IBSEB\] cost: (.*)", open(log).read())
    ok &= report("cost line reported", len(tl) > 0, tl[-1].strip() if tl else "no cost line in the log")
    print("building set:", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
