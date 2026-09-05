#!/usr/bin/env python3
"""Check the isolated-building day.

    python3 check_isolated.py ibseb_day.csv faces/day 0.5 295.0

Reads the per-building CSV (one row a minute) and the face dumps
faces/day.step*.rank*.csv (one a minute) and asserts:
  1. the balance residual stays below 1e-3 W/m2 all day;
  2. at night every face loses longwave and before sunrise the roof, which
     sees the whole sky, is colder than the air and than the walls;
  3. after sunrise the east wall is the first face to rise above the air
     next to it;
  4. the peaks come in the order east wall, roof, west wall, with the
     roof peaking in the early afternoon (12:00-15:30 solar time);
  5. at 13:00 the south wall is warmer than the north wall;
  6. the shortwave absorbed by the core roof over the day equals the
     clear-sky formulas (declination, equation of time, hour angle,
     zenith, Bird beam, Liu-Jordan diffuse) integrated independently in
     Python at the dump times, to 1 percent;
  7. the slab energy of every face changes over the day by the
     integrated conduction, to 3 percent of the integrated |G|.
"""
import sys, glob, re, math
import numpy as np

SIGMA = 5.670374419e-8

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
    steps = np.array([int(re.search(r"\.step(\d+)$", f).group(1)) for f in files])
    return steps, [load(f) for f in files]

def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})"); return ok

# Independent solar geometry (the phase 2 checker's formulas, extended to the
# beam and diffuse so the day can be integrated).
def solar(t_utc_s, lat, lon, tz, doy, S0, tau, kd):
    g = 2 * math.pi * (doy - 1) / 365.25
    decl = 0.006918 - 0.399912*math.cos(g) + 0.070257*math.sin(g) - 0.006758*math.cos(2*g) + 0.000907*math.sin(2*g) - 0.00248*math.cos(3*g) + 0.00031*math.sin(3*g)
    eot = 229.18 * (0.000075 + 0.001868*math.cos(g) - 0.032077*math.sin(g) - 0.014615*math.cos(2*g) - 0.040849*math.sin(2*g))
    t_solar = (t_utc_s / 3600.0) % 24.0 + tz + (lon - 15.0 * tz) / 15.0 + eot / 60.0
    ha = math.radians(15.0 * (t_solar - 12.0))
    latr = math.radians(lat)
    cz = math.sin(latr) * math.sin(decl) + math.cos(latr) * math.cos(decl) * math.cos(ha)
    g2 = 2 * math.pi * (doy - 1) / 365.0
    e0 = 1.000110 + 0.034221*math.cos(g2) + 0.001280*math.sin(g2) + 0.000719*math.cos(2*g2) + 0.000077*math.sin(2*g2)
    if cz <= 1e-3: return cz, 0.0, 0.0
    tr = tau ** (1.0 / cz)
    return cz, S0 * e0 * tr, kd * S0 * e0 * cz * (1.0 - tr)

def orientation(d):
    return {"roof": d["dir"] == 2, "east": (d["dir"] == 0) & (d["side"] == -1), "west": (d["dir"] == 0) & (d["side"] == 1),
            "north": (d["dir"] == 1) & (d["side"] == -1), "south": (d["dir"] == 1) & (d["side"] == 1)}

def main():
    csv, prefix, dt, T_int = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    c = np.genfromtxt(csv, delimiter=",", names=True)
    hours = c["time_s"] / 3600.0
    ok = report("balance residual all day", c["resid_max_Wm2"].max() < 1e-3, f"max {c['resid_max_Wm2'].max():.1e} W/m2 over {len(c)} rows")

    steps, dumps = load_steps(prefix)
    th = (steps + 1) * dt / 3600.0
    sel = orientation(dumps[0])
    T = {k: np.array([d["T_skin"][m].mean() for d in dumps]) for k, m in sel.items()}
    Ta = np.array([d["T_air"][sel["roof"]].mean() for d in dumps])
    LW = np.array([d["LW_net"].max() for d in dumps])
    print("  hour   T_air  roof   east   south  west   north  [K]")
    for h in range(0, 24, 2):
        n = np.argmin(np.abs(th - h))
        print(f"  {th[n]:5.1f}  {Ta[n]:6.2f} " + " ".join(f"{T[k][n]:6.2f}" for k in ("roof", "east", "south", "west", "north")))

    # 2. night
    night = th < 4.5
    # The roof sees the whole sky and the walls half of it, so with the
    # faces started at the air temperature the roof ends the night coldest.
    walls_0430 = np.mean([T[k][night][-1] for k in ("east", "west", "north", "south")])
    ok &= report("night: every face loses longwave; before sunrise the roof is colder than the air and than the walls",
                 LW[night].max() < 0 and T["roof"][night][-1] < Ta[night][-1] - 1.0 and T["roof"][night][-1] < walls_0430,
                 f"max LW_net at night {LW[night].max():.1f} W/m2, at 04:30 roof - air {T['roof'][night][-1] - Ta[night][-1]:.2f} K, roof - walls {T['roof'][night][-1] - walls_0430:.2f} K")
    # 3. onset: every face is colder than the air at dawn; the first to rise
    # above the air next to it is the one the beam reaches first. (A
    # threshold on the face's own warming would pick the roof, which is the
    # coldest at dawn and rebounds on diffuse light before the beam clears
    # the horizon.)
    Tair = {k: np.array([d["T_air"][m].mean() for d in dumps]) for k, m in sel.items()}
    n0 = np.where(night)[0][-1]
    onset = {k: (th[n0:][np.where(T[k][n0:] > Tair[k][n0:])[0][0]] if (T[k][n0:] > Tair[k][n0:]).any() else 99) for k in T}
    first = min(onset, key=onset.get)
    ok &= report("east wall is the first face to rise above the air after sunrise", first == "east", ", ".join(f"{k} {onset[k]:.2f} h" for k in ("east", "roof", "north", "south", "west")))
    # 4. peaks
    peak = {k: th[np.argmax(T[k])] for k in T}
    ok &= report("peaks in the order east, roof, west with the roof in the early afternoon",
                 peak["east"] < peak["roof"] < peak["west"] and 12.0 <= peak["roof"] <= 15.5,
                 ", ".join(f"{k} {peak[k]:.2f} h ({T[k].max():.1f} K)" for k in ("east", "roof", "west", "south", "north")))
    # 5. midday contrast
    n13 = np.argmin(np.abs(th - 13.0))
    ok &= report("south wall warmer than the north wall at 13:00", T["south"][n13] > T["north"][n13] + 1.0,
                 f"south {T['south'][n13]:.2f} K, north {T['north'][n13]:.2f} K")
    # 6. daily shortwave on the core roof against the independent clear sky
    d0 = dumps[0]
    core = sel["roof"] & (d0["z_m"] > d0["z_m"][sel["roof"]].max() - 1.0)
    alb = d0["albedo"][core].mean(); fs = d0["f_sky"][core]; fg = d0["f_ground"][core]
    E_erf = sum(d["SW_abs"][core].mean() for d in dumps) * dt * (steps[1] - steps[0])
    E_py = 0.0
    for s in steps:
        t = 25200.0 + (s + 1) * dt
        cz, dni, dif = solar(t, 40.0, -105.0, -7.0, 172, 1361.0, 0.7, 0.5)
        direct = dni * max(cz, 0.0)
        diffuse = (fs * dif + fg * 0.2 * (direct + dif)).mean()
        E_py += (1.0 - alb) * (direct + diffuse) * dt * (steps[1] - steps[0])
    ok &= report("daily shortwave absorbed by the core roof vs the clear-sky formulas integrated in Python",
                 abs(E_erf - E_py) < 0.01 * E_py, f"ERF {E_erf/1e6:.3f} MJ/m2, Python {E_py/1e6:.3f} MJ/m2")
    # 7. slab energy over the day
    N = 0
    while f"T_slab_{N}" in d0: N += 1
    S0 = np.stack([dumps[0][f"T_slab_{l}"] for l in range(N)], 1); S1 = np.stack([dumps[-1][f"T_slab_{l}"] for l in range(N)], 1)
    k = d0["k_therm"]; dz = d0["thickness"] / N; rc = d0["rho_cp"]
    dE = (rc * dz) * (S1 - S0).sum(1)
    interval = dt * (steps[1] - steps[0])
    Gnet = np.zeros_like(dE); Gabs = np.zeros_like(dE)
    for a, b in zip(dumps[:-1], dumps[1:]):
        Gb_a = 2 * k * (a[f"T_slab_{N-1}"] - T_int) / dz; Gb_b = 2 * k * (b[f"T_slab_{N-1}"] - T_int) / dz
        Gnet += 0.5 * interval * ((a["G"] - Gb_a) + (b["G"] - Gb_b)); Gabs += 0.5 * interval * (np.abs(a["G"]) + np.abs(b["G"]))
    err = np.abs(dE - Gnet).max() / Gabs.max()
    ok &= report("slab energy over the day vs the integrated conduction (trapezoid on minute dumps)", err < 0.03,
                 f"worst face {err*100:.2f} % of the integrated |G| ({Gabs.max()/1e6:.2f} MJ/m2)")
    print("isolated building:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
