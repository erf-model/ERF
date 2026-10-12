#!/usr/bin/env python3
"""Check the building faces that take their radiation from the two-stream columns
(erf.ibseb.radiation = two_stream).

    check_ibseb_two_stream_provider.py <work dir> --nz <cells in z> --cosz <cos zenith>
        --tau <SW optical depth per layer> --two-level-cosz <cos zenith> --two-level-ref-z <ratio>
        --toa <W/m2> --albedo <ground albedo> --ssa <single-scattering albedo>
        --tau-lw <LW optical depth per layer> --ground-temp <K> --ground-emissivity <e>
        [--precision single]

Tests/RunIBSEBTwoStreamProvider.cmake runs the legs below under <work dir>; this script
reads their face dumps (faces.stepNNNNNN.rank*.csv, faces/set[.lev1].rank*.csv) and
building reports (ibseb_buildings.csv).
  two_stream/   a 40 m cube under a transparent sky, two steps, then restarted from step 1
  prescribed/   the same cube on the faces' own clear-sky radiation set to that sky (the same
                sun, no diffuse light, no sky longwave, the columns' ground), also restarted
  absorbing/    shortwave optical depth tau per layer, no scattering, so the beam at
                interface m is S0 cos z exp(-tau (nz - m) / cos z); a longwave optical
                depth per layer and a ground warmer than the air (--ground-temp)
  scattering/   the absorbing sky with a share --ssa of its extinction scattering
  night/        the sun below the horizon, a 70 m tower beside the cube, traps on
  two_level_*/  IBSEB_TwoStreamProviderTwoLevel.i: a cube on level 1, a taller tower outside
                it, level 1 refined in z too (--two-level-ref-z times nz layers), the
                columns' longwave on the faces; a transparent and an absorbing sky

Each face reads its column at its own height: a roof at the bottom of its fluid cell
(interface m = k), a wall at the mean of its cell's two interfaces (k and k + 1).
The checks, and the defect each one catches:
  1. Every leg dumped the same faces: the cube's roof and four walls, some sunlit, some
     shaded, some that see the ground. (Nothing to compare is a failure.)
  2. Transparent sky: every face gets the same shortwave and longwave as with the faces' own
     clear-sky radiation, at both steps. Catches a beam not divided by cos z, the total
     light taken as the diffuse light, and the ground's reflection or emission missing.
  3. Absorbing sky: every sunlit face's beam is the transparent one times
     exp(-tau (nz - m) / cos z) at its own height; the light the ground reflects onto a
     wall is times exp(-tau nz / cos z) (the beam down to the ground) times exp(-2 tau m)
     (the reflected light back up to the face; the solver passes exp(-2 tau) per layer
     without scattering). Catches one height for every face, and the reflection read at the
     ground instead of at the face.
  4. Night: the roofs' sky longwave (LW_ext / f_sky; a roof sees no ground) is positive and
     lower on the tower than on the cube, as less air lies above it. Catches the longwave
     read at one height for every face, or not at all.
  5. Before the first sweep (the report of step 0) the faces have no radiation.
  6. Two levels: on each level the beam and the reflected light follow check 3 with that
     level's own layer count (nz on level 0, r nz on level 1), and every roof gets sky
     longwave. Catches a refined level reading at a height of its own or with the layer count
     of level 0. (That each face reads its own column, not another one, the unit test
     IBSEBTwoStreamFaces checks, on columns that all differ.)
  7. Night: no shortwave on any face, and every roof gets sky longwave; the run did not trap.
  8. Restart: for both kinds of radiation, the step the restart starts from keeps its one
     report row and its dump. Catches a restart that reports that step again (with no
     radiation yet, or with the sun of the step's end).
  9. Scattering: every face's beam is that of the absorbing sky (scattering takes light out
     of the beam no more than absorption does), and every face's diffuse light is
     f_sky (down - beam) + f_ground up at its height, with the column solved here
     independently (the two-stream layer solution and adding method of
     ERF_TwoStreamSW.H, written out below). Catches diffuse light leaking into the beam, the
     diffuse sky dropped, scaled or read at another height.
 10. Longwave at the face's height: under the absorbing sky, with the ground warmer than the
     air, every face's sky and ground longwave is f_sky down + f_ground up at its height in
     a gray column rebuilt here from the faces' air temperatures (linear in height for the
     deck's neutral sounding). Catches the walls reading the ground's longwave at the
     ground instead of at their height (17 W/m2 on the higher walls here).
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


# Relative tolerances of the comparisons: round-off of ERF's Real (--precision).
RTOL = 1e-9


def close(a, b, rtol=None, atol=1e-9):
    rtol = RTOL if rtol is None else rtol
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def at_sample(key, t):
    """t at a face's sample: a roof's interface k, a wall's mean of k and k + 1."""
    k, d = key[2], key[3]
    return t(k) if d == 2 else 0.5 * (t(k) + t(k + 1))


def beam_factor(key, nz, tau, mu):
    """Transmission of the beam from the top of the atmosphere down to a face's sample."""
    return at_sample(key, lambda m: math.exp(-tau * (nz - m) / mu))


def reflection_factor(key, nz, tau, mu):
    """Transmission of the light the ground reflects, as a face sees it: the beam down to the
    ground, then the upwelling up to the face's sample. With no scattering each layer passes
    exp(-2 tau) of the upwelling (the two-stream diffusivity factor 2)."""
    return math.exp(-tau * nz / mu) * at_sample(key, lambda m: math.exp(-2.0 * tau * m))



def sw_layer(tau, w0, mu0, g=0.0):
    """Reflection and transmission of one homogeneous layer (ERF_TwoStreamSW.H: the
    Eddington-type gamma coefficients, diffuse and direct-beam parts)."""
    g1 = (8 - w0 * (5 + 3 * g)) / 4
    g2 = 3 * w0 * (1 - g) / 4
    g3 = (2 - 3 * g * mu0) / 4
    g4 = 1 - g3
    a1, a2 = g1 * g4 + g2 * g3, g1 * g3 + g2 * g4
    k = math.sqrt(max((g1 - g2) * (g1 + g2), 1e-12))
    e = math.exp(-k * tau)
    e2 = e * e
    rt = 1 / (k * (1 + e2) + g1 * (1 - e2))
    r_dif, t_dif = rt * g2 * (1 - e2), rt * 2 * k * e
    kmu = k * mu0
    if abs(1 - kmu * kmu) < 1e-4:
        kmu = 1 - 1e-2 if kmu < 1 else 1 + 1e-2
    t_n = math.exp(-tau / mu0)
    rt2 = w0 * rt / (1 - kmu * kmu)
    kg3, kg4 = k * g3, k * g4
    r_dir = max(rt2 * ((1 - kmu) * (a2 + kg3) - (1 + kmu) * (a2 - kg3) * e2 - 2 * (kg3 - a2 * kmu) * e * t_n), 0)
    t_dir = max(-rt2 * ((1 + kmu) * (a1 + kg4) * t_n - (1 - kmu) * (a1 - kg4) * e2 * t_n
                        - 2 * (kg4 + a1 * kmu) * e), 0)
    if r_dir + t_dir > 1 - t_n:
        f = (1 - t_n) / (r_dir + t_dir)
        r_dir, t_dir = r_dir * f, t_dir * f
    return r_dif, t_dif, r_dir, t_dir


def sw_column(nz, tau, w0, mu, toa, alb):
    """Beam, total down and up at interfaces 0 .. nz of a uniform column over a ground of
    albedo alb, by the adding method (ERF_TwoStreamColumn.H)."""
    fdir = [toa * mu * math.exp(-tau * (nz - m) / mu) for m in range(nz + 1)]
    r_dif, t_dif, r_dir, t_dir = sw_layer(tau, w0, mu)
    a, src = [alb] + [0.0] * nz, [alb * fdir[0]] + [0.0] * nz
    for m in range(nz):
        den = max(1 - r_dif * a[m], 1e-12)
        a[m + 1] = r_dif + t_dif * t_dif * a[m] / den
        src[m + 1] = r_dir * fdir[m + 1] + t_dif * (src[m] + a[m] * t_dir * fdir[m + 1]) / den
    up, dn = [0.0] * (nz + 1), [0.0] * (nz + 1)
    up[nz], dn[nz] = src[nz], fdir[nz]
    d_above = 0.0
    for m in range(nz - 1, -1, -1):
        den = max(1 - r_dif * a[m], 1e-12)
        d = (t_dif * d_above + t_dir * fdir[m + 1] + r_dif * src[m]) / den
        up[m], dn[m], d_above = a[m] * d + src[m], fdir[m] + d, d
    return fdir, dn, up


def lw_column(temps, tau_lw, t_ground, eps_ground):
    """Gray longwave down and up at interfaces 0 .. nz (ERF_TwoStreamLW.H), no longwave
    coming in at the top, the ground emitting and reflecting."""
    sig, t = 5.670374419e-8, math.exp(-tau_lw)
    nz = len(temps)
    dn = [0.0] * (nz + 1)
    for m in range(nz - 1, -1, -1):
        dn[m] = dn[m + 1] * t + sig * temps[m] ** 4 * (1 - t)
    up = [eps_ground * sig * t_ground ** 4 + (1 - eps_ground) * dn[0]] + [0.0] * nz
    for m in range(nz):
        up[m + 1] = up[m] * t + sig * temps[m] ** 4 * (1 - t)
    return dn, up


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("--nz", type=int, required=True)
    ap.add_argument("--cosz", type=float, required=True)
    ap.add_argument("--tau", type=float, required=True)
    ap.add_argument("--two-level-cosz", type=float, required=True)
    ap.add_argument("--two-level-ref-z", type=int, required=True)
    ap.add_argument("--toa", type=float, required=True)
    ap.add_argument("--albedo", type=float, required=True)
    ap.add_argument("--ssa", type=float, required=True)
    ap.add_argument("--tau-lw", type=float, required=True)
    ap.add_argument("--ground-temp", type=float, required=True)
    ap.add_argument("--ground-emissivity", type=float, required=True)
    ap.add_argument("--precision", choices=("double", "single"), default="double")
    a = ap.parse_args()
    global RTOL
    RTOL = 1e-9 if a.precision == "double" else 2e-5
    ok = True

    def check(name, cond, detail):
        nonlocal ok
        print(f"  {name}: {'PASS' if cond else 'FAIL'} ({detail})")
        ok &= bool(cond)

    w = a.work
    ts, pre, ab, sc = (step_dump(f"{w}/{r}", 1) for r in ("two_stream", "prescribed", "absorbing", "scattering"))
    ts2, pre2 = step_dump(f"{w}/two_stream", 2), step_dump(f"{w}/prescribed", 2)
    sunlit = [k for k, f in ts.items() if f["SW_direct_in"] > 1.0]
    shaded = [k for k, f in ts.items() if f["shadow"] > 0.5]
    ground = [k for k, f in ts.items() if f["f_ground"] > 0.05]
    dirs = sorted({k[3] for k in ts})
    check("1. faces to compare", len(ts) > 0 and set(ts) == set(pre) == set(ab) == set(sc) == set(ts2) == set(pre2)
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
    check("2. transparent sky: the same as the faces' own clear-sky radiation", worst is None,
          f"every face and term within {RTOL:g} at steps 1 and 2" if worst is None else f"first mismatch {worst}")

    refl = math.exp(-a.tau * a.nz / a.cosz)
    bad_dir = [k for k in sunlit
               if not close(ab[k]["SW_direct_in"], ts[k]["SW_direct_in"] * beam_factor(k, a.nz, a.tau, a.cosz), RTOL)]
    bad_dif = [k for k in ts
               if not close(ab[k]["SW_diffuse_in"], ts[k]["SW_diffuse_in"] * reflection_factor(k, a.nz, a.tau, a.cosz), RTOL)]
    heights = sorted({k[2] for k in sunlit})
    check("3. absorbing sky: the beam at each face's height, the reflection from the ground",
          len(heights) > 2 and not bad_dir and not bad_dif,
          f"{len(sunlit)} sunlit faces in {len(heights)} cells of height, beam to the ground {refl:.6f}; "
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
        nz = a.nz * (a.two_level_ref_z if lev == 1 else 1)
        lit = [k for k, f in c.items() if f["SW_direct_in"] > 1.0]
        bad = [k for k in lit
               if not close(b[k]["SW_direct_in"], c[k]["SW_direct_in"] * beam_factor(k, nz, a.tau, mu), RTOL)]
        bad += [k for k in c
                if not close(b[k]["SW_diffuse_in"], c[k]["SW_diffuse_in"] * reflection_factor(k, nz, a.tau, mu), RTOL)]
        cold = [k for k in c if k[3] == 2 and not (c[k]["LW_ext"] > 0.0 and b[k]["LW_ext"] > 0.0)]
        good &= bool(lit) and not bad and not cold
        detail.append(f"level {lev} ({nz} layers): {len(c)} faces, {len(lit)} sunlit, {len(bad)} mismatches, "
                      f"{len(cold)} roofs without sky longwave")
    check("6. two levels, level 1 refined in z: each level's columns at each face's height", good, "; ".join(detail))

    lit_n = [k for k, f in night.items() if f["SW_abs"] != 0.0]
    cold_n = [k for k, f in night.items() if f["dir"] == 2 and not f["LW_ext"] > 0.0]
    check("7. night: no shortwave, the sky's longwave", night and not lit_n and not cold_n,
          f"{len(night)} faces at step 2, {len(lit_n)} with shortwave, {len(cold_n)} roofs without sky longwave")

    good, detail = True, []
    for leg in ("two_stream", "prescribed"):
        with open(f"{w}/{leg}/ibseb_buildings.csv") as f:
            rows = list(csv.DictReader(f))
        s1 = [r for r in rows if r["step"] == "1"]
        s2 = [r for r in rows if r["step"] == "2"]
        dump1 = step_dump(f"{w}/{leg}", 1)
        good &= (len(s1) == 1 and float(s1[0]["SW_abs_mean_Wm2"]) > 0.0 and len(s2) == 2
                 and any(f["SW_abs"] > 0.0 for f in dump1.values()))
        detail.append(f"{leg}: {len(s1)} step-1 row(s), SW_abs_mean {[r['SW_abs_mean_Wm2'] for r in s1]}, "
                      f"{len(s2)} step-2 rows (before and after the restart), step-1 dump SW_abs max "
                      f"{max((f['SW_abs'] for f in dump1.values()), default=float('nan')):.1f}")
    check("8. a restart keeps the report and dump of the step it starts from", good, "; ".join(detail))

    fdir, dn, up = sw_column(a.nz, a.tau, a.ssa, a.cosz, a.toa, a.albedo)
    beam_diff = [k for k in ts if not close(sc[k]["SW_direct_in"], ab[k]["SW_direct_in"], RTOL)]
    expect = {k: f["f_sky"] * at_sample(k, lambda m: dn[m] - fdir[m]) + f["f_ground"] * at_sample(k, lambda m: up[m])
              for k, f in sc.items()}
    bad = [k for k in sc if not close(sc[k]["SW_diffuse_in"], expect[k], RTOL)]
    roofs = [k for k in sc if k[3] == 2 and sc[k]["f_sky"] > 0.5]
    check("9. scattering sky: the same beam, the diffuse light of the column at each face's height",
          roofs and min(sc[k]["SW_diffuse_in"] for k in roofs) > 1.0 and not beam_diff and not bad,
          f"{len(beam_diff)} faces whose beam changed; {len(bad)} of {len(sc)} faces off the solved column"
          + (f" (first {bad[0]}: {sc[bad[0]]['SW_diffuse_in']} vs {expect[bad[0]]})" if bad else "")
          + f"; roof diffuse {min((sc[k]['SW_diffuse_in'] for k in roofs), default=0):.1f}-"
          f"{max((sc[k]['SW_diffuse_in'] for k in roofs), default=0):.1f} W/m2")

    # The air temperature of the column, linear in height for the neutral sounding, from the
    # walls' fluid cells; the longwave column rebuilt on it.
    walls = [(k[2], f["T_air"]) for k, f in ab.items() if k[3] != 2]
    kbar = sum(k for k, _ in walls) / len(walls)
    tbar = sum(t for _, t in walls) / len(walls)
    slope = sum((k - kbar) * (t - tbar) for k, t in walls) / sum((k - kbar) ** 2 for k, _ in walls)
    temps = [tbar + slope * (k - kbar) for k in range(a.nz)]
    ldn, lup = lw_column(temps, a.tau_lw, a.ground_temp, a.ground_emissivity)
    lw_tol = max(RTOL, 1e-6)
    lexp = {k: f["f_sky"] * at_sample(k, lambda m: ldn[m]) + f["f_ground"] * at_sample(k, lambda m: lup[m])
            for k, f in ab.items()}
    lbad = [k for k in ab if not close(ab[k]["LW_ext"], lexp[k], lw_tol)]
    at_ground = max((f["f_ground"] * (lup[0] - at_sample(k, lambda m: lup[m])) for k, f in ab.items()), default=0.0)
    check("10. longwave at each face's height, the ground warmer than the air",
          ab and not lbad and at_ground > 1.0,
          f"{len(lbad)} of {len(ab)} faces off the rebuilt column (relative {lw_tol:g})"
          + (f" (first {lbad[0]}: {ab[lbad[0]]['LW_ext']} vs {lexp[lbad[0]]})" if lbad else "")
          + f"; reading the ground's longwave at the ground would add up to {at_ground:.1f} W/m2")

    print("ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
