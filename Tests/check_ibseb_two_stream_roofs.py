#!/usr/bin/env python3
"""Checks of erf.ibseb.two_stream_roofs (CTest IBSEB_TwoStreamRoofs): the two-stream column
under a building starts on the building's roof.

With the input, a column under a building is the column of the layers above the roof, on
a surface with the roof face's albedo, emissivity and skin temperature; without it the
column passes through the building's cells down to the held ground under it. The expected
values come from the column solvers of check_ibseb_two_stream_provider.py (the adding
method and the grey longwave written again in Python) run on the shorter column, and from
the face dumps.

 1. the legs ran with the same faces, and the cube has roofs;
 2. shortwave: the diffuse light on every roof face is that of the shorter column over a
    surface of the roof's albedo with the input, and of the whole column over the ground
    without it (times the roof's sky view where a taller part of the building stands
    beside it). A roof the taller part shades reflects no beam in its column, so less of
    its own reflection comes back. The beam on the roofs is the same in both runs;
 3. longwave under a sky that neither absorbs nor emits: the outgoing longwave over every
    building column is emissivity sigma T_roof^4 of its own roof face with the input (the
    skin temperature the faces held before the sweep), and the held ground's without it;
    under an absorbing sky, that of the grey column of the layers above the roof over the
    roof's temperature and emissivity (the air's temperatures from the wall faces), and,
    as that sky hides the surface from the top, the longwave the walls above a lower roof
    take from the roof's column at their own height;
 4. nothing changes in the open columns in the first step: the outgoing longwave and the
    heating rates of every open column, and every wall face beside open ground, are those
    of the run without the input, to the last bit; the walls that rise above a lower roof
    read that roof's column and do change;
 5. the air above the roofs is heated differently, and the building's cells not at all;
 6. the ground under the building stays held;
 7. one rank against several;
 8. a restart;
 9. two levels: on the refined level the roofs' diffuse light changes and that of the walls
    beside open ground does not;
10. with erf.ibseb.ground_view_shadow too: the walls above a lower roof take the roof as
    the ground of their column (a sunlit ground view above zero, at most the ground view).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_ibseb_two_stream_buildings as bld2  # noqa: E402  (the plotfile reader)
import check_ibseb_two_stream_provider as prov  # noqa: E402  (the dump reader and the column solvers)

SIGMA = 5.670374419e-8


def fields(run, prefix, step, lev=0, two_d=False):
    names, levels, _ = bld2.read_plotfile(f"{run}/{prefix}{step:05d}")
    return {n: bld2.field(levels, lev, n, two_d=two_d) for n in names}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    for name, kind in (("nx", int), ("ny", int), ("nz", int), ("cosz", float), ("toa", float), ("tau", float),
                       ("ssa", float), ("tau-lw", float), ("albedo", float), ("emissivity", float), ("roof-albedo", float),
                       ("roof-emissivity", float)):
        ap.add_argument(f"--{name}", type=kind, required=True)
    ap.add_argument("--precision", choices=("double", "single"), default="double")
    a = ap.parse_args()
    rtol = 1e-9 if a.precision == "double" else 5e-5
    ok = True

    def check(name, cond, detail):
        nonlocal ok
        print(f"  {name}: {'PASS' if cond else 'FAIL'} ({detail})")
        ok &= bool(cond)

    def close(x, y):
        return abs(x - y) <= rtol * (abs(y) + 1.0)

    w = a.work
    legs = ("clear_off", "clear_on", "scattering_off", "scattering_on", "serial")
    dump0 = {leg: prov.step_dump(f"{w}/{leg}", 0) for leg in legs}
    dump1 = {leg: prov.step_dump(f"{w}/{leg}", 1) for leg in legs}
    ref = dump1["scattering_on"]
    same = all(set(d) == set(ref) for d in list(dump0.values()) + list(dump1.values()))
    roofs = [k for k in ref if k[3] == 2 and k[4] < 0]
    walls = [k for k in ref if k[3] != 2]
    # The cell on each building column's roof: the roof face's own cell.
    roof_k = {(k[0], k[1]): k[2] for k in roofs}
    heights = sorted(set(roof_k.values()))
    check("1. the same faces in every leg, the cube's roofs",
          ref and same and len(roofs) > 10 and len(walls) > 10 and len(heights) > 1 and 0 < heights[0] and heights[-1] < a.nz,
          f"{len(ref)} faces, the same in {len(legs)} legs: {same}; {len(roofs)} roof faces on cell {heights}, {len(walls)} wall faces")
    if not ok:
        print("FAILED")
        return 1
    cols = [(i, j) for i in range(a.nx) for j in range(a.ny)]
    bcols = sorted(roof_k)
    ocols = [c for c in cols if c not in roof_k]

    # 2. The diffuse light on the roofs: total down less the beam at the column's lower boundary.
    # One column per roof height (the cube's edge columns are lower than its core).
    fd_off, dn_off, _ = prov.sw_column(a.nz, a.tau, a.ssa, a.cosz, a.toa, a.albedo)
    exp_on, exp_off = {}, {}
    for kr in heights:
        for lit in (0.0, 1.0):
            fd_on, dn_on, _ = prov.sw_column(a.nz - kr, a.tau, a.ssa, a.cosz, a.toa, a.roof_albedo, ground_beam=lit)
            exp_on[(kr, lit)] = dn_on[0] - fd_on[0]
        exp_off[kr] = dn_off[kr] - fd_off[kr]
    on, off = dump1["scattering_on"], dump1["scattering_off"]
    # A roof below a taller part of the building sees less than the whole sky (f_sky < 1)
    # and no ground; where that part shades it (the face's own shadow), its column reflects
    # no beam.
    shaded = [k for k in roofs if on[k]["shadow"] > 0.5]
    bad_on = [k for k in roofs if not close(on[k]["SW_diffuse_in"], on[k]["f_sky"] * exp_on[(k[2], 1.0 - on[k]["shadow"])])]
    bad_off = [k for k in roofs if not close(off[k]["SW_diffuse_in"], off[k]["f_sky"] * exp_off[k[2]])]
    partial = [k for k in roofs if on[k]["f_sky"] < 0.999]
    beam_same = all(on[k]["SW_direct_in"] == off[k]["SW_direct_in"] for k in roofs) and \
        max(on[k]["SW_direct_in"] for k in roofs) > 100.0
    top = heights[-1]
    check("2. the diffuse light on the roofs",
          not bad_on and not bad_off and beam_same and len(shaded) > 2 and len(shaded) < len(roofs)
          and abs(exp_on[(top, 1.0)] - 57.4153) < 1e-3 and abs(exp_off[top] - 52.8873) < 1e-3
          and exp_on[(top, 0.0)] < exp_on[(top, 1.0)] - 5.0,
          f"{len(roofs) - len(bad_on)} of {len(roofs)} roofs follow the column above the roof over the roof's albedo with "
          f"the input ({exp_on[(top, 1.0)]:.4f} W/m2 on cell {top} in the sun, {exp_on[(top, 0.0)]:.4f} in a shadow), "
          f"{len(roofs) - len(bad_off)} the whole column over the ground without ({exp_off[top]:.4f} W/m2); roofs on cells "
          f"{heights}, {len(partial)} of them under a taller part, {len(shaded)} in its shadow; the same beam on the roofs: {beam_same}")

    # 3. The outgoing longwave under a sky without longwave: the surface's own emission.
    olr_on = fields(f"{w}/clear_on", "plt2d", 1, two_d=True)["OLR"]
    olr_off = fields(f"{w}/clear_off", "plt2d", 1, two_d=True)["OLR"]
    t_g0 = fields(f"{w}/clear_off", "plt2d", 0, two_d=True)["seb_t_sfc"]
    skin0 = {(k[0], k[1]): f["T_skin"] for k, f in dump0["clear_on"].items() if k in roofs}
    emis = {(k[0], k[1]): f["emissivity"] for k, f in dump0["clear_on"].items() if k in roofs}
    bad_on = [c for c in bcols if not close(olr_on[c], emis[c] * SIGMA * skin0[c] ** 4)]
    bad_off = [c for c in bcols if not close(olr_off[c], a.emissivity * SIGMA * t_g0[c] ** 4)]
    gap = max(abs(olr_on[c] - olr_off[c]) for c in bcols)
    check("3. outgoing longwave over the building, transparent sky",
          not bad_on and not bad_off and gap > 1.0 and all(close(emis[c], a.roof_emissivity) for c in bcols),
          f"{len(bcols) - len(bad_on)} of {len(bcols)} building columns at emissivity sigma T_roof^4 with the input "
          f"({olr_on[bcols[0]]:.3f} W/m2), {len(bcols) - len(bad_off)} at the held ground's without ({olr_off[bcols[0]]:.3f} W/m2)")

    # 3b. Under the absorbing sky: the grey column of the layers above the roof. The air's
    # temperatures are those the wall faces report (a straight line in the height here).
    air = [(k[2], f["T_air"]) for k, f in dump1["scattering_on"].items() if k[3] != 2]
    kbar, tbar = sum(k for k, _ in air) / len(air), sum(t for _, t in air) / len(air)
    slope = sum((k - kbar) * (t - tbar) for k, t in air) / sum((k - kbar) ** 2 for k, _ in air)
    temps = [tbar + slope * (k - kbar) for k in range(a.nz)]
    olr_abs = fields(f"{w}/scattering_on", "plt2d", 1, two_d=True)["OLR"]
    skin_s = {(k[0], k[1]): f["T_skin"] for k, f in dump0["scattering_on"].items() if k in roofs}
    lw_tol = max(rtol, 1e-6)
    bad_lw, sample = [], None
    for c in bcols:
        _, up = prov.lw_column(temps[roof_k[c]:], a.tau_lw, skin_s[c], emis[c])
        if abs(olr_abs[c] - up[-1]) > lw_tol * up[-1]:
            bad_lw.append(c)
        sample = up[-1]
    # This sky hides the surface from the top (the two runs differ there by millimetres of a
    # W/m2), so the roof is also checked from close by: a wall that rises above a lower roof
    # takes, through its sky and ground views, the longwave coming down and going up in
    # the roof's column at its own height.
    on_s, off_s = dump1["scattering_on"], dump1["scattering_off"]
    over = [k for k in on_s if k[3] != 2 and (k[0], k[1]) in roof_k]
    bad_wall, apart = [], []
    for k in over:
        c = (k[0], k[1])
        dn, up = prov.lw_column(temps[roof_k[c]:], a.tau_lw, skin_s[c], emis[c])
        expect = on_s[k]["f_sky"] * prov.at_sample(k, lambda m: dn[m - roof_k[c]]) \
            + on_s[k]["f_ground"] * prov.at_sample(k, lambda m: up[m - roof_k[c]])
        if abs(on_s[k]["LW_ext"] - expect) > lw_tol * (expect + 1.0):
            bad_wall.append(k)
        apart.append(abs(on_s[k]["LW_ext"] - off_s[k]["LW_ext"]))
    check("3b. longwave over the building, absorbing sky",
          # (the two runs must stand apart by several times the bound the walls are held to,
          # in either precision)
          not bad_lw and over and not bad_wall and max(apart) > 4 * lw_tol * sample,
          f"{len(bcols) - len(bad_lw)} of {len(bcols)} building columns send out the longwave of the grey column of the layers "
          f"above their roof ({sample:.3f} W/m2 on the last); {len(over) - len(bad_wall)} of the {len(over)} walls above a lower "
          f"roof take that column's longwave at their height, up to {max(apart):.2f} W/m2 from the run without the input")

    # 4-6. The scattering legs, first step.
    p_on, p_off = fields(f"{w}/scattering_on", "plt", 1), fields(f"{w}/scattering_off", "plt", 1)
    q_on, q_off = fields(f"{w}/scattering_on", "plt2d", 1, two_d=True), fields(f"{w}/scattering_off", "plt2d", 1, two_d=True)
    open_same = all(q_on["OLR"][c] == q_off["OLR"][c] for c in ocols) and \
        all(p_on[n][(c[0], c[1], k)] == p_off[n][(c[0], c[1], k)] for n in ("qsrc_sw", "qsrc_lw") for c in ocols for k in range(a.nz))
    wall_cols = ("SW_direct_in", "SW_diffuse_in", "LW_in", "SW_abs")
    walls_open = [k for k in walls if (k[0], k[1]) not in roof_k]
    walls_roof = [k for k in walls if (k[0], k[1]) in roof_k]
    walls_same = all(on[k][c] == off[k][c] for k in walls_open for c in wall_cols)
    walls_moved = [k for k in walls_roof if abs(on[k]["SW_diffuse_in"] - off[k]["SW_diffuse_in"]) > 1e-6]
    olr_moved = min(abs(q_on["OLR"][c] - q_off["OLR"][c]) for c in bcols)
    check("4. the open columns and the walls in the first step",
          open_same and walls_same and olr_moved > 1e-3 and len(ocols) > 900 and len(walls_open) > 10
          and walls_roof and len(walls_moved) == len(walls_roof),
          f"{len(ocols)} open columns with the same outgoing longwave and heating rates: {open_same}; {len(walls_open)} wall "
          f"faces beside open ground with the same light and longwave: {walls_same}; {len(walls_moved)} of the "
          f"{len(walls_roof)} walls above a lower roof change; the outgoing longwave over the building moves by at "
          f"least {olr_moved:.3f} W/m2")
    above = [(c[0], c[1], k) for c in bcols for k in range(roof_k[c], a.nz)]
    inside = [(c[0], c[1], k) for c in bcols for k in range(roof_k[c])]
    d_sw = max(abs(p_on["qsrc_sw"][x] - p_off["qsrc_sw"][x]) / abs(p_off["qsrc_sw"][x]) for x in above)
    d_lw = max(abs(p_on["qsrc_lw"][x] - p_off["qsrc_lw"][x]) / (abs(p_off["qsrc_lw"][x]) + 1e-30) for x in above)
    inside_zero = all(p[n][x] == 0.0 for p in (p_on, p_off) for n in ("qsrc_sw", "qsrc_lw") for x in inside)
    check("5. heating above the roofs and inside the building",
          0.15 < d_sw < 0.25 and 0.03 < d_lw < 0.08 and inside_zero,
          f"the shortwave heating above the roofs changes by up to {100 * d_sw:.1f} %, the longwave by up to "
          f"{100 * d_lw:.1f} %; no heating in the {len(inside)} cells of the building: {inside_zero}")
    held = all(q_on["seb_t_sfc"][c] == q_off["seb_t_sfc"][c] and q_on["seb_sw_net"][c] == 0.0 for c in bcols)
    check("6. the ground under the building", held, f"{len(bcols)} columns held, no shortwave absorbed: {held}")

    # 7. One rank against several.
    s_one = fields(f"{w}/serial", "plt2d", 1, two_d=True)
    p_one = fields(f"{w}/serial", "plt", 1)
    worst = max(max(abs(s_one[n][c] - q_on[n][c]) / (abs(q_on[n][c]) + 1.0) for n in ("OLR", "seb_t_sfc", "seb_lw_net") for c in cols),
                max(abs(p_one[n][x] - p_on[n][x]) / (abs(p_on[n][x]) + 1e-12) for n in ("qsrc_sw", "qsrc_lw") for x in above))
    check("7. one rank against several", worst <= 100 * rtol, f"largest relative difference {worst:.1e}")

    # 8. A restart from the checkpoint after the first step.
    c2, r2 = fields(f"{w}/scattering_on", "plt2d", 2, two_d=True), fields(f"{w}/restart", "plt2d", 2, two_d=True)
    pc2, pr2 = fields(f"{w}/scattering_on", "plt", 2), fields(f"{w}/restart", "plt", 2)
    worst = max(max(abs(c2[n][c] - r2[n][c]) / (abs(c2[n][c]) + 1.0) for n in ("OLR", "seb_t_sfc", "seb_lw_net") for c in cols),
                max(abs(pc2[n][x] - pr2[n][x]) / (abs(pc2[n][x]) + 1e-12) for n in ("qsrc_sw", "qsrc_lw") for x in above))
    moved = max(abs(c2["OLR"][c] - q_on["OLR"][c]) for c in bcols)
    check("8. a restart", worst <= 100 * rtol and moved > 1e-6,
          f"largest relative difference at step 2 {worst:.1e}; the outgoing longwave over the roofs moved by up to "
          f"{moved:.4f} W/m2 in that step (the roofs warm)")

    # 9. Two levels.
    f_on = prov.faces(f"{w}/two_level/faces/set.lev1.rank*.csv")
    f_off = prov.faces(f"{w}/two_level_off/faces/set.lev1.rank*.csv")
    roofs1 = [k for k in f_on if k[3] == 2 and k[4] < 0]
    walls1 = [k for k in f_on if k[3] != 2]
    roof_k1 = {(k[0], k[1]) for k in roofs1}
    walls1 = [k for k in walls1 if (k[0], k[1]) not in roof_k1]   # the walls beside open ground
    roofs_moved = [k for k in roofs1 if abs(f_on[k]["SW_diffuse_in"] - f_off[k]["SW_diffuse_in"]) > 1e-6]
    walls_same1 = all(f_on[k]["SW_diffuse_in"] == f_off[k]["SW_diffuse_in"] for k in walls1)
    check("9. two levels",
          set(f_on) == set(f_off) and len(roofs1) > 10 and len(roofs_moved) == len(roofs1) and walls_same1,
          f"{len(f_on)} faces on level 1; the diffuse light changes on {len(roofs_moved)} of {len(roofs1)} roofs and on "
          f"none of the {len(walls1)} walls beside open ground: {walls_same1}")

    # 10. With erf.ibseb.ground_view_shadow as well: a wall above a lower roof reads the
    # beam and the reflection at the roof its column starts on, not at a ground that
    # column no longer reaches.
    g = prov.step_dump(f"{w}/shadow_on", 1)
    over_roof = [k for k in g if k[3] != 2 and (k[0], k[1]) in roof_k and g[k]["f_ground"] > 0.0]
    lit_some = [k for k in over_roof if g[k]["f_ground_lit"] > 0.0]
    # The bounds follow --precision (rtol): in single precision the two view fractions are
    # sums of ray weights in float, equal only to its roundoff.
    bounded = all(-rtol <= g[k]["f_ground_lit"] <= g[k]["f_ground"] + rtol for k in g if k[3] != 2)
    less = all(g[k]["SW_diffuse_in"] <= on[k]["SW_diffuse_in"] * (1.0 + rtol) + rtol for k in g if k[3] != 2)
    # Where part of such a wall's ground view is shaded, it loses the beam's share of the
    # light going up at the roof: strictly less diffuse light than without the option.
    # A shaded roof reflects no beam, so there the share is zero and nothing is lost.
    roof_lit = {(k[0], k[1]): 1.0 - g[k]["shadow"] for k in g if k[3] == 2 and k[4] < 0}
    part = [k for k in over_roof if g[k]["f_ground_lit"] < g[k]["f_ground"] - rtol]
    part_lit = [k for k in part if roof_lit[(k[0], k[1])] > 0.5]
    part_dark = [k for k in part if roof_lit[(k[0], k[1])] < 0.5]
    loses = [k for k in part_lit if g[k]["SW_diffuse_in"] < on[k]["SW_diffuse_in"] - 1e-3]
    keeps = [k for k in part_dark if close(g[k]["SW_diffuse_in"], on[k]["SW_diffuse_in"])]
    check("10. with the shadow on the ground the walls see",
          set(g) == set(on) and over_roof and len(lit_some) > 0.5 * len(over_roof) and bounded and less
          and len(part) > 3 and len(loses) == len(part_lit) and len(keeps) == len(part_dark),
          f"{len(lit_some)} of the {len(over_roof)} walls above a lower roof that see ground have a sunlit ground view above "
          f"zero; within the ground view on every wall: {bounded}; of the {len(part)} with part of that view shaded, "
          f"{len(loses)} of the {len(part_lit)} above a sunlit roof take less diffuse light than without the option and "
          f"{len(keeps)} of the {len(part_dark)} above a shaded roof the same; none takes more: {less}")

    print("PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
