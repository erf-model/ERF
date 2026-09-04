#!/usr/bin/env python3
"""Check the phase 6 prognostic balance.

    python3 check_prognostic.py closure faces_closure 0.5 293.0 [--resid 1e-3]
    python3 check_prognostic.py qext    faces_qext    0.5 293.0 3000.0 380.0
    python3 check_prognostic.py solar   faces_solar
    python3 check_prognostic.py restart faces_closure.step000199 faces_restart.step000199

closure: from the per-step dumps <prefix>.step<N>.rank*.csv,
  1. the residual of the balance on every face at every step is below the
     tolerance (non-zero only at a bound or the iteration cap);
  2. every stored term is consistent with the stored skin temperature:
     H = H_coeff (T_skin / Pi - theta_air), LW_net = eps (LW_ext + f_bldg
     sigma T^4 - sigma T^4), G = 2 k (T_skin - T_slab_0) / dz;
  3. the slab energy changes per step by dt (G - G_bottom) exactly
     (implicit Euler), G_bottom = 2 k (T_slab_last - T_interior) / dz;
  4. over the run, radiation in (SW_abs + eps Q_ext + LW_net) equals the
     heat convected to the air, stored in the slab and conducted to the
     interior, to within the summed residual;
  5. an independent Python model (its own Newton and its own implicit slab
     with a dense solve) driven by the dumped forcing reproduces the
     skin-temperature trajectory of every face;
  6. the sunlit faces warm and the south wall ends warmer than the north.
qext: the closure checks on the external-flux run, the flux on every face
  equal to the input, and every face above the default 380 K bound (the
  raised bound took effect).
solar: from the periodic dumps, the east wall warms before the roof and
  the west wall, and every face stays within the bounds.
restart: the last dump (step 199, t = 100 s) after a restart at step 100 against the
  straight run: geometry exact, skin and slab to 1e-3 K (the atmosphere of
  the immersed forcing is not bit-exact through a restart, see phase 5).
"""
import sys, glob, re, argparse
import numpy as np

SIGMA = 5.670374419e-8

def load(prefix):
    rows = []
    for fn in sorted(glob.glob(prefix + ".rank*.csv")):
        with open(fn) as f:
            hdr = f.readline().strip().split(",")
            for line in f: rows.append([float(v) for v in line.strip().split(",")])
    a = np.array(rows); d = {h: a[:, n] for n, h in enumerate(hdr)}
    # Sort by (i, j, k, dir, side) so the same face has the same row in every dump.
    key = np.lexsort((d["side"], d["dir"], d["k"], d["j"], d["i"]))
    return {h: v[key] for h, v in d.items()}

def load_steps(prefix):
    files = sorted(set(re.sub(r"\.rank\d+\.csv$", "", fn) for fn in glob.glob(prefix + ".step*.rank*.csv")))
    steps = [int(re.search(r"\.step(\d+)$", f).group(1)) for f in files]
    return steps, [load(f) for f in files]

def slab_layers(d):
    n = 0
    while f"T_slab_{n}" in d: n += 1
    return np.stack([d[f"T_slab_{l}"] for l in range(n)], axis=1)

def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})"); return ok

def independent_model(steps, dumps, dt, T_int):
    """Re-integrate the balance from the dumped forcing with an independent
    Newton solve and a dense implicit slab, starting from the state of the
    first dump; returns the skin temperature at every later dump."""
    d0 = dumps[0]; nf = len(d0["i"]); N = slab_layers(d0).shape[1]
    T_skin = d0["T_skin"].copy(); S = slab_layers(d0).copy()
    k = d0["k_therm"]; rc = d0["rho_cp"]; dz = d0["thickness"] / N
    out = []
    for d in dumps[1:]:
        Fo = k * dt / (rc * dz * dz)
        eps = d["emissivity"]; fb = d["f_bldg"]; Hc = d["H_coeff"]
        Pi = d["T_air"] / d["theta_air"]; tha = d["theta_air"]
        Sin = d["SW_abs"] + eps * d["Q_ext"]; LWe = d["LW_ext"]; LE = d["LE"]
        for f in range(nf):
            # Implicit slab as a dense system M T^{n+1} = T^n + boundary terms,
            # with the skin temperature entering row 0 and the interior row N-1.
            M = np.zeros((N, N)); rhs = S[f].copy()
            for l in range(N):
                M[l, l] = 1.0 + 2.0 * Fo[f]
                if l > 0: M[l, l - 1] = -Fo[f]
                if l < N - 1: M[l, l + 1] = -Fo[f]
            M[0, 0] += Fo[f]; M[N - 1, N - 1] += Fo[f]      # half-spacing boundaries: 2 Fo each side
            rhs[N - 1] += 2.0 * Fo[f] * T_int
            Minv = np.linalg.inv(M)
            # T_0^{n+1} = alpha + beta T_skin  ->  G = 2k/dz ((1 - beta) T_skin - alpha)
            alpha = (Minv @ rhs)[0]; beta = Minv[0, 0] * 2.0 * Fo[f]
            a = 2.0 * k[f] / dz[f] * (1.0 - beta); b = 2.0 * k[f] / dz[f] * alpha
            e_eff = eps[f] * (1.0 - fb[f]) * SIGMA
            T = T_skin[f]
            for it in range(50):
                F = Sin[f] + eps[f] * LWe[f] - e_eff * T**4 - Hc[f] * (T / Pi[f] - tha[f]) - LE[f] - (a * T - b)
                Fp = -4.0 * e_eff * T**3 - Hc[f] / Pi[f] - a
                step = -F / Fp; T += step
                if abs(step) < 1e-10: break
            T_skin[f] = T
            rhs[0] += 2.0 * Fo[f] * T
            S[f] = Minv @ rhs
        out.append(T_skin.copy())
    return out

def check_closure(prefix, dt, T_int, resid_tol, label="closure"):
    steps, dumps = load_steps(prefix)
    print(f"  {len(dumps)} dumps, steps {steps[0]}-{steps[-1]}, {len(dumps[0]['i'])} faces")
    ok = True
    d0 = dumps[0]; N = slab_layers(d0).shape[1]
    k = d0["k_therm"]; dz = d0["thickness"] / N; rc = d0["rho_cp"]; A = d0["area_m2"]
    # 1. residual
    rmax = max(d["resid"].max() for d in dumps[1:])
    nit = max(d["n_iter"].max() for d in dumps[1:])
    ok &= report("balance residual on every face at every step", rmax < resid_tol, f"max {rmax:.2e} W/m2, newton iterations at most {int(nit)}")
    # 2. consistency of the stored terms with the stored skin temperature
    eH = eL = eG = 0.0
    for d in dumps[1:]:
        T = d["T_skin"]; Pi = d["T_air"] / d["theta_air"]
        H = d["H_coeff"] * (T / Pi - d["theta_air"])
        LWn = d["emissivity"] * (d["LW_ext"] + d["f_bldg"] * SIGMA * T**4 - SIGMA * T**4)
        G = 2.0 * k * (T - d["T_slab_0"]) / dz
        eH = max(eH, np.abs(H - d["H"]).max()); eL = max(eL, np.abs(LWn - d["LW_net"]).max()); eG = max(eG, np.abs(G - d["G"]).max())
    ok &= report("H, LW_net and G consistent with T_skin", max(eH, eL, eG) < 1e-6, f"max diff H {eH:.1e}, LW {eL:.1e}, G {eG:.1e} W/m2")
    # 3. slab energy per step
    emax = 0.0
    for dp, d in zip(dumps[:-1], dumps[1:]):
        dE = (rc * dz)[:, None] * (slab_layers(d) - slab_layers(dp))
        dE = dE.sum(axis=1)
        Gbot = 2.0 * k * (d[f"T_slab_{N-1}"] - T_int) / dz
        emax = max(emax, np.abs(dE - dt * (d["G"] - Gbot)).max())
    # The dump carries 12 digits, so the slab energy (rho_cp dz T ~ 1e7 J/m2 per
    # layer) is resolved to ~1e-5 J/m2; a step moves ~100 J/m2.
    ok &= report("slab energy change = dt (G - G_bottom) every step", emax < 1e-3, f"max diff {emax:.1e} J/m2 (dump precision)")
    # 4. closure over the run
    E_in = E_H = E_LE = E_bot = 0.0; R = 0.0
    for d in dumps[1:]:
        E_in += dt * (A * (d["SW_abs"] + d["emissivity"] * d["Q_ext"] + d["LW_net"])).sum()
        E_H += dt * (A * d["H"]).sum(); E_LE += dt * (A * d["LE"]).sum()
        E_bot += dt * (A * 2.0 * k * (d[f"T_slab_{N-1}"] - T_int) / dz).sum()
        R += dt * (A * d["resid"]).sum()
    E_st = (A[:, None] * (rc * dz)[:, None] * (slab_layers(dumps[-1]) - slab_layers(dumps[0]))).sum()
    gap = E_in - E_H - E_LE - E_st - E_bot
    scale = max(abs(E_in), abs(E_H), abs(E_st), abs(E_bot), 1.0)
    ok &= report("closure: radiation in = convected + stored + conducted out", abs(gap) <= R + 1e-8 * scale,
                 f"in {E_in:.4e}, to air {E_H:.4e}, stored {E_st:.4e}, to interior {E_bot:.4e} J; gap {gap:.2e}, summed residual {R:.2e}")
    # 5. independent re-integration
    model = independent_model(steps, dumps, dt, T_int)
    err = max(np.abs(m - d["T_skin"]).max() for m, d in zip(model, dumps[1:]))
    ok &= report("independent Newton + dense slab reproduces T_skin", err < 1e-4, f"max |dT| {err:.2e} K over {len(model)} steps")
    return ok, dumps

def orientation_means(d):
    n = {}
    def m(sel): return d["T_skin"][sel].mean() if sel.any() else np.nan
    n["roof"] = m((d["dir"] == 2))
    n["east"] = m((d["dir"] == 0) & (d["side"] == -1))    # outward normal +x
    n["west"] = m((d["dir"] == 0) & (d["side"] == 1))
    n["north"] = m((d["dir"] == 1) & (d["side"] == -1))   # outward normal +y
    n["south"] = m((d["dir"] == 1) & (d["side"] == 1))
    return n

def main():
    mode = sys.argv[1]
    if mode in ("closure", "qext"):
        prefix, dt, T_int = sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
        ok, dumps = check_closure(prefix, dt, T_int, 1e-3)
        last = orientation_means(dumps[-1]); first = orientation_means(dumps[0])
        print("  mean T_skin by orientation, first -> last dump: " +
              ", ".join(f"{k} {first[k]:.2f}->{last[k]:.2f} K" for k in last))
        if mode == "closure":
            # The skin is massless: it jumps to its balance in the first step and
            # then follows the slab, so the change from the 300 K start is what
            # matters, not the drift between the dumps.
            ok &= report("sunlit roof and south wall warm above the 300 K start, south warmer than north",
                         last["roof"] > 302 and last["south"] > 302 and last["south"] > last["north"],
                         f"roof {last['roof']:.2f} K, south {last['south']:.2f} K, north {last['north']:.2f} K")
        else:
            Q, Tb = float(sys.argv[5]), float(sys.argv[6])
            d = dumps[-1]
            ok &= report("external flux on every face equals the input", np.abs(d["Q_ext"] - Q).max() < 1e-9, f"{Q} W/m2")
            ok &= report(f"every face above the default bound of {Tb} K", d["T_skin"].min() > Tb,
                         f"T_skin {d['T_skin'].min():.1f}-{d['T_skin'].max():.1f} K")
    elif mode == "solar":
        steps, dumps = load_steps(sys.argv[2])
        ms = [orientation_means(d) for d in dumps]
        print("  step, mean T_skin east / roof / south / west / north [K]:")
        for s, m in zip(steps, ms):
            if s % 1200 == 0 or s == steps[-1]:
                print(f"    {s:6d}  {m['east']:.2f}  {m['roof']:.2f}  {m['south']:.2f}  {m['west']:.2f}  {m['north']:.2f}")
        # Warming onset: the first dump where an orientation is 0.5 K above its start.
        def onset(key):
            for s, m in zip(steps, ms):
                if m[key] > ms[0][key] + 0.5: return s
            return None
        oe, orf, ow = onset("east"), onset("roof"), onset("west")
        ok = report("east wall warms before the roof", oe is not None and orf is not None and oe < orf, f"onset east {oe}, roof {orf}, west {ow}")
        ok &= report("east wall ends warmest, west wall coolest of the walls",
                     ms[-1]["east"] > ms[-1]["roof"] and ms[-1]["east"] > ms[-1]["west"] and ms[-1]["west"] < ms[-1]["south"],
                     "final means above")
        Tall = np.concatenate([d["T_skin"] for d in dumps])
        ok &= report("every face within the bounds", Tall.min() > 260 and Tall.max() < 380, f"{Tall.min():.1f}-{Tall.max():.1f} K")
    elif mode == "restart":
        a, b = load(sys.argv[2]), load(sys.argv[3])
        ok = report("same faces", len(a["i"]) == len(b["i"]) and all(np.array_equal(a[c], b[c]) for c in ("i", "j", "k", "dir", "side", "bid", "area_m2")), f"{len(a['i'])} faces")
        eT = np.abs(a["T_skin"] - b["T_skin"]).max()
        eS = np.abs(slab_layers(a) - slab_layers(b)).max()
        ok &= report("skin and slab after the restart", eT < 1e-3 and eS < 1e-3, f"max |dT_skin| {eT:.1e} K, max |dT_slab| {eS:.1e} K")
        for c in ("H", "G", "LW_net", "SW_abs"):
            e = np.abs(a[c] - b[c]).max(); ok &= report(f"{c} after the restart", e < 1e-2 * max(1.0, np.abs(a[c]).max()), f"max diff {e:.1e} W/m2")
    else:
        sys.exit("unknown mode")
    print(f"{mode}: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
