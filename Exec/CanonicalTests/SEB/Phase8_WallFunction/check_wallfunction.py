#!/usr/bin/env python3
"""Check the wall function beyond neutral (phase 8).

    python3 check_wallfunction.py neutral   faces_neutral faces_deardorff
    python3 check_wallfunction.py deardorff faces_deardorff 0.5 1000.0 1.2
    python3 check_wallfunction.py stability faces_stability faces_deardorff 1000.0 1.2
    python3 check_wallfunction.py bulkri    run_bulkri.log faces_bulkri 95.0

neutral:   a hot roof in calm air sheds almost nothing with the neutral log
           law (the tangential wind sits at its 1e-3 m/s floor) and at least
           ten times more with the convective velocity scale.
deardorff: on every face of the last dump, w* equals the Deardorff scale
           built from the previous dump's H and the depth (the mixed layer
           above a roof, z_i - z_face floored at the building height; the
           building height for a wall), and u* and H follow the neutral
           log law on sqrt(U_tan^2 + (beta w*)^2), to 1e-6.
stability: on the roofs the Obukhov length is negative and consistent with
           the face's u* and theta* (from the previous step's skin, which
           is what the wall function saw), and u* and H equal the log law
           with Dyer's psi_m and psi_h at delta / L (the code's functions),
           to 1e-6; the roof flux exceeds the run without the functions.
bulkri:    the depth diagnosed on the capped sounding is the first cell
           centre above the inversion, and the roofs' depth in w* is that
           minus the roof height.
"""
import sys, glob, re, math
import numpy as np

KAPPA = 0.41; G = 9.81; CP = 1004.5

def load(prefix):
    rows = []
    for fn in sorted(glob.glob(prefix + ".rank*.csv")):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            for line in f: rows.append([float(v) for v in line.strip().split(",")])
    a = np.array(rows); d = {h: a[:, n] for n, h in enumerate(hdr)}
    key = np.lexsort((d["side"], d["dir"], d["k"], d["j"], d["i"]))
    return {h: v[key] for h, v in d.items()}

def last_two(prefix):
    files = sorted(set(re.sub(r"\.rank\d+\.csv$", "", fn) for fn in glob.glob(prefix + ".step*.rank*.csv")))
    return load(files[-2]), load(files[-1])

def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})"); return ok

def psi_m(z):
    if z > 0: return -5.0 * z
    x = math.sqrt(math.sqrt(1.0 - 16.0 * z))
    return 2 * math.log(0.5 * (1 + x)) + math.log(0.5 * (1 + x * x)) - 2 * math.atan(x) + math.pi / 2
def psi_h(z):
    if z > 0: return -5.0 * z
    x = math.sqrt(1.0 - 16.0 * z); return 2 * math.log(0.5 * (1 + x))

def log_law(d, Ueff, lnm_corr=None, lnh_corr=None):
    """u* and H from the neutral (or corrected) log law on the dumped state."""
    delta = np.where(d["dir"] == 2, 5.0, 5.0)     # half a 10 m cell
    lnm = np.log(delta / 0.01); lnh = np.log(delta / 0.001)
    if lnm_corr is not None: lnm = lnm - lnm_corr
    if lnh_corr is not None: lnh = lnh - lnh_corr
    ustar = KAPPA * Ueff / lnm
    th_skin = d["T_skin"] * d["theta_air"] / d["T_air"]
    thstar = KAPPA * (d["theta_air"] - th_skin) / lnh
    return ustar, -d["rho"] * CP * ustar * thstar

def main():
    mode = sys.argv[1]; ok = True
    if mode == "neutral":
        a = last_two(sys.argv[2])[1]; b = last_two(sys.argv[3])[1]
        ra = a["dir"] == 2; rb = b["dir"] == 2
        Ha, Hb = a["H"][ra].mean(), b["H"][rb].mean()
        ok &= report("roof in calm air: neutral law sheds little, the convective scale ten times more",
                     Ha < 5.0 and Hb > 10 * max(Ha, 0.1), f"roof H neutral {Ha:.2f} W/m2 at T_skin {a['T_skin'][ra].mean():.1f} K, with w* {Hb:.1f} W/m2 at {b['T_skin'][rb].mean():.1f} K; U_tan {a['U_tan'][ra].max():.1e} m/s")
    elif mode == "deardorff":
        prev, d = last_two(sys.argv[2]); dt = float(sys.argv[3]); z_i = float(sys.argv[4]); beta = float(sys.argv[5])
        depth = np.where(d["dir"] == 2, np.maximum(z_i - d["z_m"], d["h_bld"]), d["h_bld"])
        Hp = np.maximum(prev["H"], 0.0)
        wstar = np.cbrt(G / d["theta_air"] * Hp / (d["rho"] * CP) * depth)
        ok &= report("depth in w*: mixed layer above the roofs, building height on the walls", np.abs(d["z_i"] - depth).max() < 1e-9,
                     f"roofs {d['z_i'][d['dir']==2].min():.0f} m, walls {d['z_i'][d['dir']<2].max():.0f} m")
        e = np.abs(d["w_star"] - wstar).max()
        ok &= report("w* is the Deardorff scale of the previous step's flux", e < 1e-9, f"max diff {e:.1e} m/s, w* up to {d['w_star'].max():.2f} m/s")
        Ueff = np.sqrt(d["U_tan"] ** 2 + (beta * d["w_star"]) ** 2)
        us, H = log_law(d, Ueff)
        e1 = np.abs(us - d["ustar"]).max() / d["ustar"].max(); e2 = np.abs(H - d["H"]).max() / np.abs(d["H"]).max()
        ok &= report("u* and H from the log law on sqrt(U_tan^2 + (beta w*)^2)", e1 < 1e-6 and e2 < 1e-6, f"rel diff u* {e1:.1e}, H {e2:.1e}")
    elif mode == "stability":
        prev, d = last_two(sys.argv[2]); b = last_two(sys.argv[3])[1]; z_i = float(sys.argv[4]); beta = float(sys.argv[5])
        roof = d["dir"] == 2
        L = d["Olen"][roof]
        ok &= report("roof Obukhov length negative (unstable)", (L < 0).all() and (np.abs(L) < 1e5).all(), f"L in [{L.min():.1f}, {L.max():.1f}] m")
        # The wall function runs before the balance, on the skin temperature of
        # the previous step; the stored H is the balance's flux at the new one.
        # So theta* for L comes from the previous dump's skin, with the
        # corrected heat log law at the stored L.
        us = d["ustar"][roof]; th = d["theta_air"][roof]
        th_skin_old = prev["T_skin"][roof] * th / d["T_air"][roof]
        lnh_eff = np.log(5.0 / 0.001) - np.array([psi_h(5.0 / l) for l in L])
        thstar = KAPPA * (th - th_skin_old) / lnh_eff
        Lc = us * us * th / (KAPPA * G * thstar)
        e = np.abs(Lc - L).max() / np.abs(L).max()
        ok &= report("L = u*^2 theta / (kappa g theta*) with the face's own u* and theta* (previous skin)", e < 1e-6, f"rel diff {e:.1e}")
        Ueff = np.sqrt(d["U_tan"] ** 2 + (beta * d["w_star"]) ** 2)
        zeta = 5.0 / d["Olen"]
        pm = np.array([psi_m(z) if abs(z) < 1e4 else 0.0 for z in zeta]); ph = np.array([psi_h(z) if abs(z) < 1e4 else 0.0 for z in zeta])
        pm[d["dir"] < 2] = 0.0; ph[d["dir"] < 2] = 0.0      # walls stay neutral
        usl, Hl = log_law(d, Ueff, pm, ph)
        e1 = np.abs(usl - d["ustar"]).max() / d["ustar"].max(); e2 = np.abs(Hl - d["H"]).max() / np.abs(d["H"]).max()
        ok &= report("u* and H from the log law with Dyer's psi_m, psi_h at delta/L on the roofs", e1 < 1e-6 and e2 < 1e-6, f"rel diff u* {e1:.1e}, H {e2:.1e}")
        ok &= report("roof flux above the run without the functions", d["H"][roof].mean() > b["H"][b["dir"] == 2].mean(),
                     f"with {d['H'][roof].mean():.1f} W/m2, without {b['H'][b['dir']==2].mean():.1f} W/m2")
    elif mode == "bulkri":
        log = open(sys.argv[2]).read(); d = last_two(sys.argv[3])[1]; z_exp = float(sys.argv[4])
        zi = [float(m) for m in re.findall(r"mixed-layer depth for w\*: ([-\d.eE+]+) m", log)]
        ok &= report("bulk Richardson depth on the capped sounding", abs(zi[-1] - z_exp) < 10.0, f"diagnosed {zi[-1]:.1f} m (expected {z_exp:.0f}), range over the run {min(zi):.1f}-{max(zi):.1f} m")
        roof = d["dir"] == 2
        ok &= report("roofs' depth in w* is the mixed layer above the roof", np.abs(d["z_i"][roof] - np.maximum(zi[-1] - d["z_m"][roof], d["h_bld"][roof])).max() < 1e-9,
                     f"{d['z_i'][roof].min():.1f} m above the {d['z_m'][roof].min():.0f} m roof")
    print(f"{mode}: {'PASS' if ok else 'FAIL'}"); sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
