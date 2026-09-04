#!/usr/bin/env python3
"""Check the phase 5 slab conduction and materials.

    python3 check_ground.py thick run_thick.log 1.5 2.0e6 20 0.25 faces_thick 0.2 200 50.0
    python3 check_ground.py thin  faces_thin 1.0 0.02 20
    python3 check_ground.py materials faces_materials materials.csv
    python3 check_ground.py restart faces_thin faces_restart

thick: the slab starts at the skin temperature (320 K) and the interior
    boundary is 20 K colder, so a thermal wave enters from the interior side.
    Over 50 s it travels a few millimetres into a 200 mm slab of 1 mm layers:
    the flux through the skin must stay zero, the top layer must stay at the
    skin temperature, and the bottom layer must follow the semi-infinite
    erfc solution for a step at a boundary.
thin: a light 20 mm slab reaches the steady linear profile within 100 s:
    G = k dT / L and the layer temperatures on the straight line.
materials: each building carries the properties of its material from the CSV.
restart: the slab and fluxes after a checkpoint restart equal the straight run.
"""
import sys, glob, math
import numpy as np

def load(prefix):
    rows = []
    for fn in sorted(glob.glob(prefix + ".rank*.csv")):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            for line in f: rows.append([float(v) for v in line.strip().split(",")])
    a = np.array(rows); return {h: a[:, n] for n, h in enumerate(hdr)}

def check_thick_log(log, k, rho_cp, dT, dt):
    import re
    t, G = [], []
    for line in open(log):
        if "[IBSEB] lev=0" in line:
            t.append(float(re.search(r" t=([-\d.eE+]+)", line).group(1)))
            G.append(float(re.search(r"G_mean=([-\d.eE+]+)", line).group(1)))
    t = np.array(t); G = np.array(G)
    alpha = k / rho_cp
    # The slab starts at the skin temperature; the interior boundary steps to
    # skin - dT at t = 0. Heat flows out of the bottom; the flux through the
    # skin stays zero until the wave crosses the slab. With the thick slab the
    # wave (sqrt(alpha t) ~ 6 mm at 50 s) never reaches the skin, so G must
    # stay near zero, and the bottom flux follows the semi-infinite solution.
    ok1 = np.all(np.abs(G[1:]) < 1.0e-3 * k * dT / 0.01)
    print(f"  skin flux stays zero while the wave from the interior is far away: {'PASS' if ok1 else 'FAIL'} (max |G| {np.abs(G[1:]).max():.2e} W/m2)")
    return ok1

def check_thick_bottom(prefix, k, rho_cp, dT, L, n, t_end):
    d = load(prefix)
    alpha = k / rho_cp; dz = L / n
    # Bottom-layer centre sits dz/2 above the interior boundary. Semi-infinite
    # solution for a step dT at the boundary: T(x, t) = T_skin - dT erfc(x / (2 sqrt(alpha t))).
    from math import erfc, sqrt
    x = dz / 2
    T_exp = 320.0 - dT * erfc(x / (2 * sqrt(alpha * t_end)))
    T_got = d["T_slab_bottom"]
    err = np.abs(T_got - T_exp).max()
    ok = err < 0.05 * dT
    print(f"  bottom layer vs semi-infinite erfc solution at t = {t_end:.1f} s: expected {T_exp:.3f} K, got {T_got.min():.3f}-{T_got.max():.3f} K -> {'PASS' if ok else 'FAIL'}")
    # And the skin-side layer is still at the skin temperature.
    ok2 = np.abs(d["T_slab_top"] - 320.0).max() < 1e-6
    print(f"  top layer still at the skin temperature: {'PASS' if ok2 else 'FAIL'}")
    return ok and ok2

def check_thin(prefix, k, L, dT):
    d = load(prefix)
    G_exp = k * dT / L
    err = np.abs(d["G"] - G_exp).max() / G_exp
    ok = err < 0.01
    print(f"  thin slab at steady state: G = k dT / L = {G_exp:.1f} W/m2, got {d['G'].min():.1f}-{d['G'].max():.1f} -> {'PASS' if ok else 'FAIL'}")
    # Linear profile: top and bottom layer centres at L/8 and 7L/8 from the skin
    n = 4; T_top = 320.0 - dT * (0.5 / n); T_bot = 320.0 - dT * ((n - 0.5) / n)
    ok2 = np.abs(d["T_slab_top"] - T_top).max() < 0.02 and np.abs(d["T_slab_bottom"] - T_bot).max() < 0.02
    print(f"  linear profile: top layer {T_top:.2f} K, bottom {T_bot:.2f} K -> {'PASS' if ok2 else 'FAIL'}")
    return ok and ok2

def check_materials(prefix, csv):
    d = load(prefix)
    mats = {}
    for line in open(csv):
        if line.startswith("mat_id") or not line.strip(): continue
        p = line.strip().split(","); mats[int(p[0])] = [float(v) for v in p[2:7]]
    ok = True
    for b, mid in ((1, 1), (2, 2)):
        m = d["bid"] == b
        vals = [d["albedo"][m], d["emissivity"][m], d["k_therm"][m], d["rho_cp"][m], d["thickness"][m]]
        good = int(d["mat"][m].min()) == mid and int(d["mat"][m].max()) == mid and all(np.allclose(v, e) for v, e in zip(vals, mats[mid]))
        ok = ok and good
        print(f"  building {b} carries material {mid} (albedo {mats[mid][0]}, k {mats[mid][2]} W/m/K, {mats[mid][4]} m): {'PASS' if good else 'FAIL'}")
    # The absorbed shortwave is zero (night) but the emissivity enters LW_net.
    return ok

def check_restart(a_prefix, b_prefix):
    """The slab-owned columns (skin, layers, G, materials, geometry) must be
    exact after the restart. The atmosphere-derived columns (wind, density,
    air temperature, H) are only required to be close: the immersed-forcing
    atmosphere of the development branch does not restart bit-for-bit (the
    wind at the faces differs by about one part in ten thousand), which is
    outside the balance."""
    def load_sorted(prefix):
        fns = sorted(glob.glob(prefix + ".rank*.csv"))
        hdr = open(fns[0]).readline().strip().split(",")
        a = np.vstack([np.loadtxt(f, delimiter=",", skiprows=1) for f in fns])
        return hdr, a[np.lexsort((a[:, 4], a[:, 3], a[:, 2], a[:, 1], a[:, 0]))]
    hdr, a = load_sorted(a_prefix); _, b = load_sorted(b_prefix)
    atm = {"T_air", "theta_air", "rho", "U_tan", "ustar", "H", "LW_in", "LW_net"}
    ok = a.shape == b.shape
    worst_exact, worst_atm = 0.0, 0.0
    for n, name in enumerate(hdr):
        d = np.abs(a[:, n] - b[:, n]).max()
        scale = 1.0 + np.abs(a[:, n]).max()
        if name in atm: worst_atm = max(worst_atm, d / scale)
        else: worst_exact = max(worst_exact, d)
    ok = ok and worst_exact < 1e-9 and worst_atm < 1e-3
    print(f"  slab, G and geometry after restart equal the straight run: {'PASS' if worst_exact < 1e-9 else 'FAIL'} (max |diff| {worst_exact:.1e})")
    print(f"  atmosphere-derived columns close (development IF restart is not bit-exact): {'PASS' if worst_atm < 1e-3 else 'FAIL'} (max rel diff {worst_atm:.1e})")
    return ok

if __name__ == "__main__":
    m = sys.argv[1]
    if m == "thick":
        ok = check_thick_log(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]))
        ok = check_thick_bottom(sys.argv[7], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[8]), int(sys.argv[9]), float(sys.argv[10])) and ok
    elif m == "thin": ok = check_thin(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
    elif m == "materials": ok = check_materials(sys.argv[2], sys.argv[3])
    else: ok = check_restart(sys.argv[2], sys.argv[3])
    sys.exit(0 if ok else 1)
