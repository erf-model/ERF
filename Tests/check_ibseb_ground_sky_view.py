#!/usr/bin/env python3
"""Checks of erf.ibseb.ground_sky_view (CTest IBSEB_GroundSkyView): the two-stream ground
balance beside the buildings.

With the input, the open ground of a column keeps the sky's diffuse light and longwave
through the part of its view the walls leave open and takes the walls' reflected light
and longwave through the rest. The walls' share of the view, v, and what they send come
from the faces' own ground rays, by the reciprocity of view factors: a face of area A_f
whose n rays of N land in a column of area A_c adds A_f n / (N A_c) to v, and that weight
times its radiosity to what the ground receives.

Everything here is recomputed from the runs' face dumps and the buildings' geometry: the
rays of every wall face are cast again as straight lines against the solid columns taken
as boxes (not the code's column walk), and the fields of the 2D plotfiles are compared
column by column.

 1. the legs ran with the same faces; the walls' share is in the plotfile only with the input;
 2. the walls' share of every column against the count from the geometry, and no ray
    on a line between two columns;
 3. the shortwave of the ground balance, input on less off, after the first step:
    (1 - alb_dif) SW_walls - v SW_dif, transparent and scattering skies;
 4. the longwave likewise: emis (LW_walls - v LW_down), with LW_down zero under the
    transparent sky and, under the scattering one, the value the run without the input
    gives from its own ground: net / emis + sigma T_ground^4;
 5. what the walls send is what the ground gets: the ground's longwave gain summed over
    the columns against the faces' radiosity times their area and ground view;
 6. the ground under the building is untouched, the shadow is the same field, and the
    ground far from the walls changes by less than 1 %;
 7. one rank against several;
 8. a restart;
 9. two levels: the refined level's share stays within what its own faces send (part of
    it lands on ground outside the level), and its ground's longwave changes exactly
    where the share is above zero;
10. a block in the periodic corner, with ray maps that reach across the edge without
    spanning the domain: the share from the geometry, and several ranks against one;
11. the walls' skin temperature is that of the start of the step: with the faces'
    temperature moving, the ground's longwave of the second step follows the skin of the
    dump after the first step, not the one after the second.
"""
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_ibseb_ground_view_shadow as gvs  # noqa: E402  (the rays and the buildings as boxes)
import check_ibseb_two_stream_buildings as bld2  # noqa: E402  (the plotfile reader)
import check_ibseb_two_stream_provider as prov  # noqa: E402  (the dump reader)

SIGMA = 5.670374419e-8


def plot2d(run, step, lev=0):
    """{name: {(i, j): value}} of one level of a 2D plotfile."""
    names, levels, _ = bld2.read_plotfile(f"{run}/plt2d{step:05d}")
    return {n: bld2.field(levels, lev, n, two_d=True) for n in names}


def walls_to_ground(faces, bld, n_az, n_el, t_skin=None):
    """Per column: the walls' share of the view, their reflected shortwave and their
    longwave arriving [W/m2]; and the number of rays on a line between two columns."""
    view, sw, lw, ties = {}, {}, {}, 0
    n_ray = n_az * n_el
    a_col = bld.dx * bld.dx
    for key, f in faces.items():
        d, side = int(f["dir"]), int(f["side"])
        if d == 2:
            continue
        p = (f["x_m"], f["y_m"], f["z_m"])
        ts = f["T_skin"] if t_skin is None else t_skin[key]["T_skin"]
        j_sw = f["albedo"] * (f["SW_direct_in"] + f["SW_diffuse_in"])
        # The longwave on the face before its balance of the step: from the sky and the
        # ground (LW_ext) and from the walls it sees, at the skin temperature it starts with.
        lw_in = f["LW_ext"] + f["f_bldg"] * SIGMA * ts ** 4
        j_lw = f["emissivity"] * SIGMA * ts ** 4 + (1.0 - f["emissivity"]) * lw_in
        w = f["area_m2"] / (n_ray * a_col)
        for ie in range(n_el):
            for ia in range(n_az):
                s = gvs.hemisphere_direction(d, -side, ia, ie, n_az, n_el)
                if s[2] >= 0.0:
                    continue
                t_g = p[2] / -s[2]
                if bld.blocked(p, s, t_g):
                    continue
                xl, yl = (p[0] + t_g * s[0]) / bld.dx, (p[1] + t_g * s[1]) / bld.dx
                if math.floor(xl - 1e-7) != math.floor(xl + 1e-7) or math.floor(yl - 1e-7) != math.floor(yl + 1e-7):
                    ties += 1
                c = (math.floor(xl) % bld.nx, math.floor(yl) % bld.ny)
                view[c] = view.get(c, 0.0) + w
                sw[c] = sw.get(c, 0.0) + w * j_sw
                lw[c] = lw.get(c, 0.0) + w * j_lw
    return view, sw, lw, ties


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    for name, kind in (("nx", int), ("ny", int), ("dx", float), ("n-az", int), ("n-el", int), ("albedo", float),
                       ("albedo-diffuse", float), ("emissivity", float)):
        ap.add_argument(f"--{name}", type=kind, required=True)
    ap.add_argument("--precision", choices=("double", "single"), default="double")
    a = ap.parse_args()
    rtol = 1e-9 if a.precision == "double" else 5e-5
    ok = True

    def check(name, cond, detail):
        nonlocal ok
        print(f"  {name}: {'PASS' if cond else 'FAIL'} ({detail})")
        ok &= bool(cond)

    def close(x, y, scale):
        return abs(x - y) <= rtol * (abs(scale) + 1.0)

    w = a.work
    legs = ("clear_off", "clear_on", "serial", "scattering_off", "scattering_on")
    dump0 = {leg: prov.step_dump(f"{w}/{leg}", 0) for leg in legs}
    dump1 = {leg: prov.step_dump(f"{w}/{leg}", 1) for leg in legs}
    plt1 = {leg: plot2d(f"{w}/{leg}", 1) for leg in legs}
    ref = dump1["clear_on"]
    same = all(set(d) == set(ref) for d in list(dump0.values()) + list(dump1.values()))
    cols = [(i, j) for i in range(a.nx) for j in range(a.ny)]
    absent = all(all(v == -999.0 for v in plt1[leg]["ibseb_ground_wall_view"].values()) for leg in ("clear_off", "scattering_off"))
    present = all(min(plt1[leg]["ibseb_ground_wall_view"].values()) >= 0.0 for leg in ("clear_on", "serial", "scattering_on"))
    check("1. the same faces in every leg, the walls' share only with the input",
          ref and same and absent and present,
          f"{len(ref)} faces, the same in {len(legs)} legs: {same}; the share -999 without the input: {absent}, set with it: {present}")
    if not ok:
        print("FAILED")
        return 1

    tops = {(k[0], k[1]): f["z_m"] for k, f in ref.items() if k[3] == 2 and k[4] < 0}
    bld = gvs.Buildings(tops, a.nx, a.ny, a.dx)
    open_cols = [c for c in cols if c not in tops]
    a_col = a.dx * a.dx

    # 2. The walls' share, from the geometry.
    view, _, _, ties = walls_to_ground(ref, bld, a.n_az, a.n_el)
    got = plt1["clear_on"]["ibseb_ground_wall_view"]
    worst = max(abs(got[c] - view.get(c, 0.0)) for c in cols)
    seen = [c for c in cols if view.get(c, 0.0) > 0.0]
    check("2. the walls' share of every column, counted from the geometry",
          worst <= (1e-9 if a.precision == "double" else 1e-5) and ties == 0 and len(seen) > 100
          and 0.25 < max(view.values()) < 0.6 and all(got[c] == 0.0 for c in tops),
          f"largest difference {worst:.2e} over {len(cols)} columns, {len(seen)} of them beside walls, largest share "
          f"{max(view.values()):.3f}; {ties} rays on a line between columns; zero under the building")

    # 3-5. The ground balance after the first step: both runs start from the same ground, so
    # the sweep's fluxes are the same and the difference is the walls' terms alone. The faces'
    # skin temperature at the start of the step is that of the dump of step 0.
    for label, off, on, alb_dif, lw_sky_zero in (("transparent sky", "clear_off", "clear_on", a.albedo, True),
                                                 ("scattering sky", "scattering_off", "scattering_on", a.albedo_diffuse, False)):
        v, sw_w, lw_w, _ = walls_to_ground(dump1[on], bld, a.n_az, a.n_el, t_skin=dump0[on])
        sw_off, sw_on = plt1[off]["seb_sw_net"], plt1[on]["seb_sw_net"]
        lw_off, lw_on = plt1[off]["seb_lw_net"], plt1[on]["seb_lw_net"]
        shadow = plt1[off]["ibseb_ground_shadow"]
        shaded = [c for c in open_cols if shadow[c] > 0.5]
        sunlit = [c for c in open_cols if shadow[c] < 0.5]
        # The absorbed diffuse light of the sweep: what the shaded ground keeps without the input.
        dif = [sw_off[c] for c in shaded]
        dif_abs = sum(dif) / len(dif)
        uniform = max(dif) - min(dif) <= rtol * (dif_abs + 1.0)
        bad_sw = [c for c in open_cols
                  if not close(sw_on[c] - sw_off[c], (1.0 - alb_dif) * sw_w.get(c, 0.0) - v.get(c, 0.0) * dif_abs, sw_off[c])]
        d_sw = [sw_on[c] - sw_off[c] for c in open_cols]
        check(f"3. shortwave of the ground, {label}",
              uniform and not bad_sw and max(d_sw) > 5.0 and (lw_sky_zero or min(d_sw) < -1.0) and len(shaded) > 5 and len(sunlit) > 100,
              f"{len(open_cols) - len(bad_sw)} of {len(open_cols)} open columns follow (1 - alb_dif) SW_walls - v SW_dif; "
              f"absorbed diffuse light {dif_abs:.2f} W/m2; change from {min(d_sw):.2f} to {max(d_sw):.2f} W/m2")
        d_lw = {c: lw_on[c] - lw_off[c] for c in open_cols}
        beside = [c for c in open_cols if v.get(c, 0.0) > 0.0]
        # LW_down from each column beside walls: emis (LW_walls - v LW_down) = d_lw.
        lw_down = [(lw_w[c] - d_lw[c] / a.emissivity) / v[c] for c in beside]
        mean_down = sum(lw_down) / len(lw_down)
        spread = max(lw_down) - min(lw_down)
        scale = max(lw_w.values())
        # The sky's longwave from the run without the input: its ground's net longwave is
        # emis (LW_down - sigma T^4), with the ground temperature the step started from.
        t0 = plot2d(f"{w}/{off}", 0)["seb_t_sfc"]
        far_col = min(open_cols, key=lambda c: v.get(c, 0.0))
        lw_down_off = lw_off[far_col] / a.emissivity + SIGMA * t0[far_col] ** 4
        if lw_sky_zero:
            good = abs(lw_down_off) <= 1e-6 and all(close(d_lw[c], a.emissivity * lw_w.get(c, 0.0), scale) for c in open_cols)
            detail = f"emis LW_walls in every open column, up to {max(d_lw.values()):.2f} W/m2 (no longwave from this sky)"
        else:
            good = spread <= 1e3 * rtol * (mean_down + 1.0) and mean_down > 50.0 and \
                abs(mean_down - lw_down_off) <= 1e3 * rtol * (lw_down_off + 1.0) and \
                all(close(d_lw[c], 0.0, 1.0) for c in open_cols if c not in beside)
            detail = (f"the sky longwave of the run without the input, {lw_down_off:.4f} W/m2, in all {len(beside)} columns "
                      f"beside walls ({mean_down:.4f} W/m2, spread {spread:.1e}); change from {min(d_lw.values()):.2f} "
                      f"to {max(d_lw.values()):.2f} W/m2")
        check(f"4. longwave of the ground, {label}", good and len(beside) > 100, detail)
        if lw_sky_zero:
            got_sum = sum(d_lw[c] for c in open_cols) * a_col
            sent = 0.0
            for key, f in dump1[on].items():
                if int(f["dir"]) == 2:
                    continue
                ts = dump0[on][key]["T_skin"]
                lw_in = f["LW_ext"] + f["f_bldg"] * SIGMA * ts ** 4
                sent += f["area_m2"] * f["f_ground"] * (f["emissivity"] * SIGMA * ts ** 4 + (1.0 - f["emissivity"]) * lw_in)
            check("5. what the walls send is what the ground gets",
                  abs(got_sum - a.emissivity * sent) <= 1e3 * rtol * sent and sent > 1e5,
                  f"the ground's longwave gain over all columns {got_sum:.6e} W, emis times the faces' radiosity over "
                  f"their ground view {a.emissivity * sent:.6e} W")
        foot_same = all(sw_on[c] == 0.0 and sw_off[c] == 0.0 and lw_on[c] == lw_off[c] for c in tops)
        shadow_same = all(plt1[on]["ibseb_ground_shadow"][c] == shadow[c] for c in cols)
        far = [c for c in sunlit if v.get(c, 0.0) < 0.005]
        far_small = all(abs(sw_on[c] - sw_off[c]) <= 0.01 * sw_off[c] for c in far)
        check(f"6. the ground under the building, the shadow, the ground far from the walls, {label}",
              foot_same and shadow_same and far_small and len(far) > 50,
              f"{len(tops)} columns under the building untouched: {foot_same}; the same shadow: {shadow_same}; "
              f"{len(far)} columns with a share below 0.005 change by less than 1 %: {far_small}")

    # 7. One rank against several.
    one_rank = plt1["serial"]
    worst = max(abs(one_rank[n][c] - plt1["clear_on"][n][c]) / (abs(one_rank[n][c]) + 1.0)
                for n in ("seb_sw_net", "seb_lw_net", "ibseb_ground_wall_view", "seb_t_sfc") for c in cols)
    check("7. one rank against several", worst <= 10 * rtol, f"largest relative difference {worst:.1e}")

    # 8. A restart from the checkpoint after the first step.
    cont, rest = plot2d(f"{w}/clear_on", 2), plot2d(f"{w}/restart", 2)
    worst = max(abs(cont[n][c] - rest[n][c]) / (abs(cont[n][c]) + 1.0)
                for n in ("seb_sw_net", "seb_lw_net", "ibseb_ground_wall_view", "seb_t_sfc") for c in cols)
    moved = max(abs(cont["seb_t_sfc"][c] - plt1["clear_on"]["seb_t_sfc"][c]) for c in cols)
    check("8. a restart", worst <= 10 * rtol and moved > 1e-6,
          f"largest relative difference at step 2 {worst:.1e}; the ground moved by up to {moved:.3f} K in that step")

    # 9. Two levels: the refined level's share against its own faces' ground view.
    two = f"{w}/two_level"
    names, levels, dxs = bld2.read_plotfile(f"{two}/plt2d00001")
    fine = bld2.field(levels, 1, "ibseb_ground_wall_view", two_d=True)
    coarse = bld2.field(levels, 0, "ibseb_ground_wall_view", two_d=True)
    faces1 = prov.faces(f"{two}/faces/set.lev1.rank*.csv")
    a_fine = dxs[1][0] * dxs[1][1]
    got_sum = sum(fine.values()) * a_fine
    sent = sum(f["area_m2"] * f["f_ground"] for f in faces1.values() if int(f["dir"]) != 2)
    _, levels_off, _ = bld2.read_plotfile(f"{w}/two_level_off/plt2d00001")
    lw1_on = bld2.field(levels, 1, "seb_lw_net", two_d=True)
    lw1_off = bld2.field(levels_off, 1, "seb_lw_net", two_d=True)
    changed = {c for c in fine if abs(lw1_on[c] - lw1_off[c]) > 1e-6}
    with_share = {c for c in fine if fine[c] > 0.0}
    check("9. two levels",
          faces1 and min(fine.values()) >= 0.0 and min(coarse.values()) >= 0.0 and max(fine.values()) > 0.1
          and max(coarse.values()) > 0.1 and 0.3 * sent <= got_sum <= sent * (1.0 + 1e3 * rtol)
          and changed == with_share and len(changed) > 100,
          f"{len(faces1)} faces on level 1; its share summed over its columns {got_sum:.1f} m2 of the faces' ground "
          f"view {sent:.1f} m2 (the rest lands on ground outside the refined level, which level 0's own faces serve); "
          f"largest share {max(fine.values()):.3f} on level 1, {max(coarse.values()):.3f} on level 0; the ground's longwave "
          f"on level 1 changes in {len(changed)} columns, those with a share above zero: {changed == with_share}")

    # 10. A block at the periodic edge of a 64 x 64 domain: on several ranks each rank's ray
    # map holds the periodic images beside its boxes; on one rank it is the whole domain.
    seam, seam1 = plot2d(f"{w}/seam", 1), plot2d(f"{w}/seam_serial", 1)
    faces_s = prov.step_dump(f"{w}/seam", 1)
    tops_s = {(k[0], k[1]): f["z_m"] for k, f in faces_s.items() if k[3] == 2 and k[4] < 0}
    bld_s = gvs.Buildings(tops_s, 64, 64, a.dx)
    view_s, _, _, ties_s = walls_to_ground(faces_s, bld_s, a.n_az, a.n_el)
    cols_s = [(i, j) for i in range(64) for j in range(64)]
    worst_geo = max(abs(seam["ibseb_ground_wall_view"][c] - view_s.get(c, 0.0)) for c in cols_s)
    worst_rank = max(abs(seam[n][c] - seam1[n][c]) / (abs(seam1[n][c]) + 1.0)
                     for n in ("seb_sw_net", "seb_lw_net", "ibseb_ground_wall_view") for c in cols_s)
    across = [c for c in cols_s if (c[0] >= 50 or c[1] >= 50) and view_s.get(c, 0.0) > 0.0]
    with open(f"{w}/seam/simulation.log") as f:
        log = f.read()
    with open(f"{w}/seam_serial/simulation.log") as f:
        log1 = f.read()
    partial = "the patches do not span the domain" in log and "every rank's patch spans the domain" in log1
    at_edge = min(c[0] for c in tops_s) == 0 and min(c[1] for c in tops_s) == 0
    check("10. a block at the periodic edge",
          worst_geo <= (1e-9 if a.precision == "double" else 1e-5) and worst_rank <= 10 * rtol and ties_s == 0
          and len(across) > 10 and partial and at_edge,
          f"largest difference from the geometry {worst_geo:.2e}, several ranks against one {worst_rank:.1e}; "
          f"{len(across)} columns across the edge see the block; ray maps short of the domain on several ranks and "
          f"spanning it on one: {partial}; the block touches the edge in x and in y: {at_edge}")

    # 11. The second step of the transparent leg, the faces' skin temperature moving. No
    # longwave comes from this sky, so the ground's net longwave is emis (LW_walls - sigma
    # T_ground^4) with the ground temperature after the first step.
    dump2 = prov.step_dump(f"{w}/clear_on", 2)
    plt2 = plot2d(f"{w}/clear_on", 2)
    t_g = plt1["clear_on"]["seb_t_sfc"]
    _, _, lw_start, _ = walls_to_ground(dump2, bld, a.n_az, a.n_el, t_skin=dump1["clear_on"])
    _, _, lw_end, _ = walls_to_ground(dump2, bld, a.n_az, a.n_el)
    walls_only = {c: plt2["seb_lw_net"][c] / a.emissivity + SIGMA * t_g[c] ** 4 for c in open_cols}
    err_start = max(abs(walls_only[c] - lw_start.get(c, 0.0)) for c in open_cols)
    err_end = max(abs(walls_only[c] - lw_end.get(c, 0.0)) for c in open_cols)
    moved_skin = max(abs(dump2[k]["T_skin"] - dump1["clear_on"][k]["T_skin"]) for k in dump2 if k[3] != 2)
    check("11. the walls' skin temperature at the start of the step",
          err_start <= 1e3 * rtol * max(lw_start.values()) and err_end > 1e3 * err_start + 1e-3 and moved_skin > 0.01,
          f"the ground's longwave from the walls in step 2 follows the skin after step 1 to {err_start:.1e} W/m2 and "
          f"differs from the skin after step 2 by up to {err_end:.3f} W/m2; the walls' skin moved by up to {moved_skin:.3f} K")

    print("PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
