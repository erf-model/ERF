#!/usr/bin/env python3
"""Check the phase 2 shortwave against an independent ray cast and analytic values.

    python3 check_shortwave.py faces_zen60 60 270 800 100 0.3 0.2
    python3 check_shortwave.py faces_solar solar <log>

Reads the per-rank face dumps <prefix>.rank*.csv. The embedded-boundary
reader turns each nodal building edge into a one-cell step, so every box is
a full-height core with a half-height rim and low corner stubs; the checker
therefore does not assume plane walls but

  1. rebuilds the column tops from the roof faces of the dump and casts a
     ray from every face toward the sun in Python (a 2D walk over columns,
     periodic like the deck), and requires the code's shadow flag to match
     on every face;
  2. requires direct = DNI max(0, n.s) (1 - shadow), the diffuse
     f_sky D + f_ground a_g (DNI cos z + D) with the view fractions read from
     the dump, and absorbed = (1 - a) x sum on every face;
  3. spot-checks the clean planes analytically: the tall core roof and the
     tall west wall are unshadowed, and the shadow on the short box's core
     west wall stops at H_core - gap tan(elevation) where the gap is from
     the tall core's east edge.
"""
import sys, glob, math, re
import numpy as np

DX = 5.0; L = 640.0; NX = 128

def load(prefix):
    rows = []
    for fn in sorted(glob.glob(prefix + ".rank*.csv")):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            for line in f:
                rows.append([float(v) for v in line.strip().split(",")])
    a = np.array(rows)
    return {h: a[:, n] for n, h in enumerate(hdr)}

def close(a, b, tol=1e-6):
    return bool(np.all(np.abs(a - b) <= tol * (1.0 + np.abs(b))))

def column_tops(d):
    top = np.zeros((NX, NX))
    roof = (d["dir"] == 2) & (d["side"] == -1)
    for i, j, z in zip(d["i"][roof].astype(int), d["j"][roof].astype(int), d["z_m"][roof]):
        top[i, j] = max(top[i, j], z)
    return top

def ray_blocked(x0, y0, z0, s, top):
    """Same model as the code, written independently: walk the columns the ray crosses."""
    sx, sy, sz = s
    if sz <= 0: return True
    zmax = top.max()
    if z0 >= zmax: return False
    t0 = 1e-3 * DX
    x, y, z = x0 + t0 * sx, y0 + t0 * sy, z0 + t0 * sz
    ci, cj = int(math.floor(x / DX)), int(math.floor(y / DX))
    stepx = 1 if sx > 0 else (-1 if sx < 0 else 0); stepy = 1 if sy > 0 else (-1 if sy < 0 else 0)
    big = 1e30
    tdx = DX / abs(sx) if stepx else big; tdy = DX / abs(sy) if stepy else big
    tmx = (((ci + 1) * DX - x) / sx) if stepx > 0 else ((ci * DX - x) / sx if stepx < 0 else big)
    tmy = (((cj + 1) * DX - y) / sy) if stepy > 0 else ((cj * DX - y) / sy if stepy < 0 else big)
    t = 0.0
    for _ in range(4 * 2 * NX + 8):
        ci %= NX; cj %= NX
        if z + 1e-6 * DX < top[ci, cj]: return True
        if tmx < tmy: t = tmx; tmx += tdx; ci += stepx
        else:         t = tmy; tmy += tdy; cj += stepy
        z = z0 + t0 * sz + t * sz
        if z >= zmax or t > 8 * L: return False
    return False

def check_fixed(prefix, zen_deg, az_deg, dni, dif, alb, alb_g):
    d = load(prefix)
    z = math.radians(zen_deg); az = math.radians(az_deg)
    s = (math.sin(z) * math.sin(az), math.sin(z) * math.cos(az), math.cos(z))
    cz, sz = math.cos(z), math.sin(z); dir_h = dni * cz
    ok = True
    def report(name, cond, extra=""):
        nonlocal ok
        print(f"  {name:52s} {'PASS' if cond else 'FAIL'} {extra}")
        ok = ok and bool(cond)
    n = np.zeros((len(d["i"]), 3))
    for f in range(len(d["i"])): n[f, int(d["dir"][f])] = -d["side"][f]
    cosi = n @ np.array(s)
    top = column_tops(d)
    # 1. shadow flag against the Python ray cast, face by face
    exp_shadow = np.zeros(len(d["i"]))
    for f in range(len(d["i"])):
        if cosi[f] > 0 and dni > 0:
            exp_shadow[f] = 1.0 if ray_blocked(d["x_m"][f], d["y_m"][f], d["z_m"][f], s, top) else 0.0
    mism = int(np.sum(exp_shadow != d["shadow"]))
    report("shadow flag = independent ray cast, every face", mism == 0, f"mismatches {mism}/{len(exp_shadow)}, shadowed {int(d['shadow'].sum())}")
    # 2. fluxes on every face
    report("direct = DNI max(0, n.s) (1 - shadow)", close(d["SW_direct_in"], dni * np.maximum(0, cosi) * (1 - d["shadow"])))
    f_sky = d["f_sky"]; f_gnd = d["f_ground"]   # sampled view fractions (phase 3), read from the dump
    report("diffuse = f_sky D + f_ground a_g (DNI cos z + D)", close(d["SW_diffuse_in"], f_sky * dif + f_gnd * alb_g * (dir_h + dif)))
    report("absorbed = (1 - albedo) (direct + diffuse)", close(d["SW_abs"], (1 - alb) * (d["SW_direct_in"] + d["SW_diffuse_in"])))
    # 3. analytic spot checks on the clean planes
    bids = np.unique(d["bid"]).astype(int)
    tops = {b: d["z_m"][d["bid"] == b].max() for b in bids}
    tall = max(tops, key=tops.get); short = min(tops, key=tops.get)
    isb = lambda b: d["bid"] == b
    west = (d["dir"] == 0) & (d["side"] == 1); east = (d["dir"] == 0) & (d["side"] == -1)
    roof = (d["dir"] == 2) & (d["side"] == -1)
    m = roof & isb(tall) & (d["z_m"] == tops[tall])
    report("tall core roof: unshadowed, direct = DNI cos z", d["shadow"][m].max() == 0 and close(d["SW_direct_in"][m], dir_h), f"n={m.sum()}")
    m = west & isb(tall)
    report("tall west wall: unshadowed, direct = DNI sin z", d["shadow"][m].max() == 0 and close(d["SW_direct_in"][m], dni * sz), f"n={m.sum()}")
    # Short core west wall: the plane with the largest x among the short box's west faces
    mw = west & isb(short); x_core = d["x_m"][mw].max(); m = mw & (d["x_m"] == x_core)
    me = east & isb(tall) & (d["z_m"] == d["z_m"][east & isb(tall)].max()); x_edge = d["x_m"][me].max()
    H_core = tops[tall]  # roof faces sit at the core top
    gap = x_core - x_edge
    h_shadow = H_core - gap * math.tan(math.pi / 2 - z)
    exp = (d["z_m"][m] < h_shadow).astype(float)
    report(f"short core west wall: shadow where z < {h_shadow:.1f} m (gap {gap:.0f} m)", bool(np.all(d["shadow"][m] == exp)), f"shadowed {int(d['shadow'][m].sum())}/{int(m.sum())}")
    return ok

def check_solar(prefix, log):
    line = [l for l in open(log) if "[IBSEB DEBUG] lev=0 sun:" in l][-1]
    zen = float(re.search(r"zenith=([-\d.eE+]+)", line).group(1))
    az  = float(re.search(r"azimuth=([-\d.eE+]+)", line).group(1))
    g = 2 * math.pi * (172 - 1) / 365.25
    decl = math.degrees(0.006918 - 0.399912*math.cos(g) + 0.070257*math.sin(g) - 0.006758*math.cos(2*g) + 0.000907*math.sin(2*g) - 0.00248*math.cos(3*g) + 0.00031*math.sin(3*g))
    zen_exp = 40.0 - decl
    ok = abs(zen - zen_exp) < 1.0 and abs(az - 180.0) < 3.0
    print(f"  solar noon Boulder: zenith {zen:.2f} (expected {zen_exp:.2f}), azimuth {az:.1f} (expected 180) -> {'PASS' if ok else 'FAIL'}")
    d = load(prefix)
    roof = d["dir"] == 2
    ok2 = d["SW_direct_in"][roof].max() > 500.0
    print(f"  clear-sky beam on roofs: max direct {d['SW_direct_in'][roof].max():.0f} W/m2 (> 500) -> {'PASS' if ok2 else 'FAIL'}")
    return ok and ok2

if __name__ == "__main__":
    if sys.argv[2] == "solar":
        ok = check_solar(sys.argv[1], sys.argv[3])
    else:
        ok = check_fixed(sys.argv[1], *[float(v) for v in sys.argv[2:8]])
    sys.exit(0 if ok else 1)
