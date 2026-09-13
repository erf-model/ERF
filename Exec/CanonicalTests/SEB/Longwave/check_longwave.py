#!/usr/bin/env python3
"""Check the phase 3 view fractions and longwave.

    python3 check_longwave.py faces_fixed fixed 300 0.9 0.95 300 16 8
    python3 check_longwave.py faces_gray  gray  0.83 0.9 0.95 300 16 8

Reads the per-rank face dumps. Checks:
  1. the three view fractions of every face sum to one;
  2. they equal an independent Python hemisphere sampling (same stratified
     directions, an independent ray walk that rebuilds the column tops from
     the roof faces) on every face;
  3. analytic values: the tall core roof sees only sky (f_sky = 1); every
     roof has f_ground = 0; the tall west wall, facing open ground, has
     f_sky = 1/2 and sees its own rim ledge a lot just above it and little
     from high up; the short box's core west wall sees the tall box;
  4. LW_in = f_sky LW_sky + f_ground eps_g sigma T_g^4 + f_bldg sigma T_s^4
     and LW_net = eps (LW_in - sigma T_s^4) on every face, with LW_sky the
     given value (fixed) or eps_sky sigma T_air^4 (gray).
"""
import sys, glob, math
import numpy as np

DX = 5.0; L = 640.0; NX = 128; SIGMA = 5.670374419e-8

def load(prefix):
    rows = []
    for fn in sorted(glob.glob(prefix + ".rank*.csv")):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            for line in f: rows.append([float(v) for v in line.strip().split(",")])
    a = np.array(rows); return {h: a[:, n] for n, h in enumerate(hdr)}

def close(a, b, tol=1e-6): return bool(np.all(np.abs(a - b) <= tol * (1.0 + np.abs(b))))

def column_tops(d):
    top = np.zeros((NX, NX)); roof = (d["dir"] == 2) & (d["side"] == -1)
    for i, j, z in zip(d["i"][roof].astype(int), d["j"][roof].astype(int), d["z_m"][roof]): top[i, j] = max(top[i, j], z)
    return top

SKY, GROUND, BUILDING = 0, 1, 2
def ray_hit(x0, y0, z0, s, top, zg=0.0):
    sx, sy, sz = s; zmax = top.max(); eps = 1e-6 * DX
    if sz > 0 and z0 >= zmax: return SKY
    t0 = 1e-3 * DX; x, y, z = x0 + t0 * sx, y0 + t0 * sy, z0 + t0 * sz
    ci, cj = int(math.floor(x / DX)), int(math.floor(y / DX))
    stepx = 1 if sx > 0 else (-1 if sx < 0 else 0); stepy = 1 if sy > 0 else (-1 if sy < 0 else 0)
    big = 1e30; tdx = DX / abs(sx) if stepx else big; tdy = DX / abs(sy) if stepy else big
    tmx = (((ci + 1) * DX - x) / sx) if stepx > 0 else ((ci * DX - x) / sx if stepx < 0 else big)
    tmy = (((cj + 1) * DX - y) / sy) if stepy > 0 else ((cj * DX - y) / sy if stepy < 0 else big)
    if stepx == 0 and stepy == 0: return SKY if sz > 0 else GROUND
    t = 0.0
    for _ in range(8 * NX + 8):
        ci %= NX; cj %= NX
        tp = top[ci, cj]; solid = tp > zg + eps
        tn = min(tmx, tmy); zout = z0 + t0 * sz + tn * sz
        if sz >= 0:
            if solid and z + eps < tp: return BUILDING
        else:
            if solid and z + eps < tp: return BUILDING
            if solid and zout < tp - eps: return BUILDING
            if zout <= zg + eps: return GROUND
        if tmx < tmy: t = tmx; tmx += tdx; ci += stepx
        else:         t = tmy; tmy += tdy; cj += stepy
        z = zout
        if sz > 0 and z >= zmax: return SKY
        if t > 8 * L: return SKY if sz >= 0 else GROUND
    return SKY if sz >= 0 else GROUND

def hemisphere(dir_, nsign, ia, ie, n_az, n_el):
    u = (ie + 0.5) / n_el; th = math.asin(math.sqrt(u)); ph = 2 * math.pi * (ia + 0.5) / n_az
    cn, ct, st = math.cos(th), math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph)
    n = [0, 0, 0]; t1 = [0, 0, 0]; t2 = [0, 0, 0]; n[dir_] = nsign
    if dir_ == 0: t1[1] = 1; t2[2] = 1
    elif dir_ == 1: t1[0] = 1; t2[2] = 1
    else: t1[0] = 1; t2[1] = 1
    return tuple(cn * n[a] + ct * t1[a] + st * t2[a] for a in range(3))

def main():
    prefix, mode = sys.argv[1], sys.argv[2]
    if mode == "fixed": lw_down, eps, eps_g, Tg, n_az, n_el = [float(v) for v in sys.argv[3:9]]
    else: eps_sky, eps, eps_g, Tg, n_az, n_el = [float(v) for v in sys.argv[3:9]]
    n_az, n_el = int(n_az), int(n_el)
    d = load(prefix); nf = len(d["i"]); ok = True
    def report(name, cond, extra=""):
        nonlocal ok; print(f"  {name:56s} {'PASS' if cond else 'FAIL'} {extra}"); ok = ok and bool(cond)
    fsum = d["f_sky"] + d["f_ground"] + d["f_bldg"]
    report("f_sky + f_ground + f_bldg = 1 on every face", close(fsum, np.ones(nf), 1e-12))
    top = column_tops(d)
    exp = np.zeros((nf, 3))
    for f in range(nf):
        c = [0, 0, 0]
        for ie in range(n_el):
            for ia in range(n_az):
                c[ray_hit(d["x_m"][f], d["y_m"][f], d["z_m"][f], hemisphere(int(d["dir"][f]), -int(d["side"][f]), ia, ie, n_az, n_el), top)] += 1
        exp[f] = np.array(c) / (n_az * n_el)
    got = np.stack([d["f_sky"], d["f_ground"], d["f_bldg"]], axis=1)
    mism = int(np.sum(np.any(np.abs(exp - got) > 1e-9, axis=1)))
    report("fractions = independent hemisphere sampling, every face", mism == 0, f"mismatches {mism}/{nf}")
    bids = np.unique(d["bid"]).astype(int); tops = {b: d["z_m"][d["bid"] == b].max() for b in bids}
    tall = max(tops, key=tops.get); short = min(tops, key=tops.get); isb = lambda b: d["bid"] == b
    roof = (d["dir"] == 2) & (d["side"] == -1); west = (d["dir"] == 0) & (d["side"] == 1)
    m = roof & isb(tall) & (d["z_m"] == tops[tall])
    report("tall core roof: f_sky = 1", close(d["f_sky"][m], np.ones(m.sum()), 1e-12), f"n={m.sum()}")
    report("every roof: f_ground = 0", np.all(d["f_ground"][roof] == 0))
    # The tall west wall faces open ground: the upward half of the hemisphere
    # is all sky. The downward half sees the rim ledge 5 m wide right below
    # the core wall, a lot from the rows just above it and little from high up.
    m = west & isb(tall) & (d["z_m"] > 70.0)
    fs = d["f_sky"][m]
    report("tall west wall above the rim: f_sky = 1/2", fs.min() > 0.45 and fs.max() < 0.55, f"f_sky {fs.min():.3f}-{fs.max():.3f}")
    # High up, the few percent that remain are rays that wrap round the
    # periodic domain and meet the box's own east wall from behind.
    hi = m & (d["z_m"] > 120.0); lo = m & (d["z_m"] < 77.0)
    report("tall west wall: ledge in view fades with height above it", d["f_bldg"][hi].max() < 0.10 and d["f_bldg"][lo].min() > 0.15,
           f"f_bldg {d['f_bldg'][hi].max():.3f} at z > 120 m, {d['f_bldg'][lo].min():.3f} just above the ledge")
    mw = west & isb(short); m = mw & (d["x_m"] == d["x_m"][mw].max())
    report("short core west wall sees the tall box: f_bldg > 0.2", d["f_bldg"][m].min() > 0.2, f"f_bldg {d['f_bldg'][m].min():.3f}-{d['f_bldg'][m].max():.3f}")
    lw_sky = np.full(nf, lw_down) if mode == "fixed" else eps_sky * SIGMA * d["T_air"] ** 4
    lw_in = d["f_sky"] * lw_sky + d["f_ground"] * eps_g * SIGMA * Tg ** 4 + d["f_bldg"] * SIGMA * d["T_skin"] ** 4
    report("LW_in = f_sky LW_sky + f_ground eps_g sigma Tg^4 + f_bldg sigma Ts^4", close(d["LW_in"], lw_in))
    report("LW_net = eps (LW_in - sigma Ts^4)", close(d["LW_net"], eps * (lw_in - SIGMA * d["T_skin"] ** 4)))
    if mode == "gray":
        report("gray sky: T_air read from the fluid cells is plausible (250-320 K)", d["T_air"].min() > 250 and d["T_air"].max() < 320, f"{d['T_air'].min():.1f}-{d['T_air'].max():.1f} K")
    sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
