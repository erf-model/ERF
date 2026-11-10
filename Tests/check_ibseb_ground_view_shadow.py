#!/usr/bin/env python3
"""Check the direct beam the ground reflects onto the building faces through the sunlit
part of each face's ground view (erf.ibseb.ground_view_shadow).

    check_ibseb_ground_view_shadow.py <work dir> --nx <columns> --ny <columns> --dx <m>
        --nz <cells in z> --n-az <view azimuths> --n-el <view elevations>
        --zenith-deg <deg> --azimuth-deg <deg> --low-zenith-deg <deg> --low-azimuth-deg <deg>
        --dni <W/m2> --diffuse <W/m2> --toa <W/m2> --albedo <ground, direct beam>
        --albedo-diffuse <ground, diffuse light> --tau <SW optical depth per layer>
        --ssa <single-scattering albedo> [--precision single]

Tests/RunIBSEBGroundViewShadow.cmake runs the legs below under <work dir>, each on the
40 m cube of Tests/test_files/IBSEB_TwoStreamProvider; this script reads their face dumps
(faces.stepNNNNNN.rank*.csv).
  prescribed_off/  the faces' own clear-sky radiation with a fixed beam and diffuse light,
                   the input off
  prescribed_on/   the same with the input on, two steps, a checkpoint after the first
  restart/         prescribed_on restarted from that checkpoint to step 2
  serial/          prescribed_on on one rank
  low_sun/         prescribed_on with the sun 18 degrees above the horizon in the west-south-west
  two_stream_on/   the two-stream columns under a transparent sky, the input on
  absorbing_off/, absorbing_on/    an absorbing sky with no scattering
  scattering_off/, scattering_on/  that sky with scattering, and a ground whose albedo for
                   diffuse light differs from its albedo for the beam
  two_level_off/, two_level_on/    IBSEB_TwoStreamProviderTwoLevel.i under the absorbing sky:
                   a cube on level 1, a taller tower outside it, the sun low in the east
                   (dumps faces/set[.lev1].rank*.csv)

The expected sunlit part of a face's ground view is counted here from the geometry, not
by the code's column walk: the solid columns are boxes from the ground to the roof faces
of the dump; each of the face's n_az x n_el hemisphere rays (the stratified, cosine-weighted
sampling of the view fractions) is intersected with those boxes and their periodic images
and with the ground; a ray that reaches the ground lands in a column, and the ray from the
centre of that column toward the sun is intersected with the boxes again. A ray that lands
within 1e-6 m of the line between two columns may be counted in either.

The checks, and the defect each one catches:
  1. Every leg dumped the same faces; the legs with the input on have the f_ground_lit
     column and the others do not (the dump of a run without the input is unchanged).
  2. Geometry: every face's ground view fraction and its sunlit part are those counted
     here, for the deck's sun and for the low sun, whose shadows are three times longer.
     Catches a shadow taken at the face instead of at the ground it sees, a landing in the
     wrong column, the wrong sun, a lost periodic image. Not vacuous: some walls see only
     sunlit ground, some both, and the low sun darkens more.
  3. The faces' own radiation: the reflected and sky light on every face is
     f_sky D + a (f_ground D + f_ground_lit B) with the input on and f_sky D +
     f_ground a (B + D) with it off (B the beam on a horizontal surface, D the diffuse
     light, a the ground's albedo). Catches the diffuse reflection shaded too, or the beam
     not shaded.
  4. Nothing else moves: the direct beam, the shadow, the view fractions and the longwave
     of every face are the same with the input on and off, and so is every roof.
  5. Two-stream columns, transparent sky: the reflected light is a S0 cos z f_ground_lit.
  6. Absorbing sky: all the light coming up is the reflected beam, so the reflected light
     with the input on is that with it off times f_ground_lit / f_ground, at every height.
  7. Scattering sky: the reflected and sky light on every face is
     f_sky (down - beam) + up (f_ground - (f_ground - f_ground_lit) s) at its height, with
     s = (direct albedo x beam at the ground) / (up at the ground), the column solved here
     (check_ibseb_two_stream_provider.py). Catches the diffuse albedo used for the beam,
     the share taken at the face instead of at the ground, the whole reflection shaded.
     Also printed: how far the light coming up is from a reference, the column solved
     again over a ground the beam does not reach, taken through the shaded part of the
     view (the sky's diffuse light kept as it is): the largest difference on a face,
     which must stay below --share-gap W/m2, and the range of the difference per unit of
     shaded view over the faces' heights.
  8. One rank and several give the same sunlit parts and the same light, to the last digit
     written. (On this small periodic domain every rank's map of the buildings is the
     whole domain; that the answer does not depend on where a map ends is the unit test
     IBSEBGroundViewShadow.SameAnswerOnAPatchAndOnTheWholeDomain.)
  9. A restart gives the dump of the step after it as the run without a restart does. The
     sunlit parts are not in the checkpoint; they are found again from the geometry.
 10. Two levels: on each level the faces have their sunlit parts, never above the ground
     view fraction; the reflected light with the input is that without it times
     f_ground_lit / f_ground (check 6, on each level's own columns); the roofs do not
     change; and the tower outside the refined level shades ground that level 1's faces
     see. Catches a refined level without its landings or with the level below's map.
Stdlib only.
"""
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_ibseb_two_stream_provider as prov  # noqa: E402  (the dump reader and the column solver)


def hemisphere_direction(d, nsign, ia, ie, n_az, n_el):
    """The (ia, ie) ray of the cosine-weighted hemisphere around the outward normal of a
    face of direction d: azimuth 2 pi (ia + 1/2) / n_az, angle from the normal
    asin(sqrt((ie + 1/2) / n_el)), tangents along the other two axes in order."""
    theta = math.asin(math.sqrt((ie + 0.5) / n_el))
    phi = 2.0 * math.pi * (ia + 0.5) / n_az
    cn, ct, st = math.cos(theta), math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi)
    n, t1, t2 = [0.0] * 3, [0.0] * 3, [0.0] * 3
    n[d] = float(nsign)
    if d == 0:
        t1[1], t2[2] = 1.0, 1.0
    elif d == 1:
        t1[0], t2[2] = 1.0, 1.0
    else:
        t1[0], t2[1] = 1.0, 1.0
    return tuple(cn * n[c] + ct * t1[c] + st * t2[c] for c in range(3))


class Buildings:
    """The solid columns as boxes from the ground to their tops, on a periodic nx x ny grid."""

    def __init__(self, tops, nx, ny, dx):
        self.tops, self.nx, self.ny, self.dx = tops, nx, ny, dx
        self.zmax = max(tops.values())

    def _entry(self, p, s, x0, y0, top):
        """Distance at which the ray from p along s enters the box of a column, or None."""
        lo, hi = (x0, y0, 0.0), (x0 + self.dx, y0 + self.dx, top)
        t_in, t_out = 0.0, float("inf")
        for c in range(3):
            if s[c] != 0.0:
                a, b = (lo[c] - p[c]) / s[c], (hi[c] - p[c]) / s[c]
                t_in, t_out = max(t_in, min(a, b)), min(t_out, max(a, b))
            elif not lo[c] < p[c] < hi[c]:
                return None
        # A face lies on the side of its own building: leaving it at once is no hit.
        return t_in if t_out > t_in + 1e-9 else None

    def blocked(self, p, s, t_end):
        """Whether the ray from p along s meets a solid column or one of its periodic images
        before the distance t_end."""
        lx, ly = self.nx * self.dx, self.ny * self.dx
        ex, ey = p[0] + t_end * s[0], p[1] + t_end * s[1]
        for (i, j), top in self.tops.items():
            x0, y0 = i * self.dx, j * self.dx
            # Only the images the ray's footprint can reach.
            ox0 = math.floor((min(p[0], ex) - x0 - self.dx) / lx)
            ox1 = math.ceil((max(p[0], ex) - x0) / lx)
            oy0 = math.floor((min(p[1], ey) - y0 - self.dx) / ly)
            oy1 = math.ceil((max(p[1], ey) - y0) / ly)
            for ox in range(ox0, ox1 + 1):
                for oy in range(oy0, oy1 + 1):
                    t = self._entry(p, s, x0 + ox * lx, y0 + oy * ly, top)
                    if t is not None and t < t_end:
                        return True
        return False

    def shaded(self, ci, cj, sun):
        """Whether the centre of the ground of column (ci, cj) is in the buildings' shadow."""
        if sun[2] <= 0.0:
            return True
        if (ci % self.nx, cj % self.ny) in self.tops:
            return True
        c = ((ci + 0.5) * self.dx, (cj + 0.5) * self.dx, 0.0)
        return self.blocked(c, sun, self.zmax / sun[2])


def ground_view(face, bld, sun, n_az, n_el, shade_cache):
    """(ground rays, sunlit ground rays at least, at most) of a face, counted from the
    geometry; the two sunlit counts differ only by rays landing on a line between columns."""
    d, side = int(face["dir"]), int(face["side"])
    p = (face["x_m"], face["y_m"], face["z_m"])
    n_gnd = lo = hi = 0
    for ie in range(n_el):
        for ia in range(n_az):
            s = hemisphere_direction(d, -side, ia, ie, n_az, n_el)
            if s[2] >= 0.0:
                continue
            t_g = p[2] / -s[2]
            if bld.blocked(p, s, t_g):
                continue
            n_gnd += 1
            xl, yl = (p[0] + t_g * s[0]) / bld.dx, (p[1] + t_g * s[1]) / bld.dx
            cols = set()
            for ci in {math.floor(xl - 1e-7), math.floor(xl + 1e-7)}:
                for cj in {math.floor(yl - 1e-7), math.floor(yl + 1e-7)}:
                    cols.add((ci, cj))
            lit = []
            for c in cols:
                if c not in shade_cache:
                    shade_cache[c] = bld.shaded(c[0], c[1], sun)
                lit.append(not shade_cache[c])
            lo += all(lit)
            hi += any(lit)
    return n_gnd, lo, hi


def sun_vector(zenith_deg, azimuth_deg):
    z, a = math.radians(zenith_deg), math.radians(azimuth_deg)
    return (math.sin(a) * math.sin(z), math.cos(a) * math.sin(z), math.cos(z))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    for name, kind in (("nx", int), ("ny", int), ("dx", float), ("nz", int), ("n-az", int), ("n-el", int),
                       ("zenith-deg", float), ("azimuth-deg", float), ("low-zenith-deg", float),
                       ("low-azimuth-deg", float), ("dni", float), ("diffuse", float), ("toa", float),
                       ("albedo", float), ("albedo-diffuse", float), ("tau", float), ("ssa", float)):
        ap.add_argument(f"--{name}", type=kind, required=True)
    ap.add_argument("--share-gap", type=float, default=0.6)
    ap.add_argument("--precision", choices=("double", "single"), default="double")
    a = ap.parse_args()
    rtol = 1e-9 if a.precision == "double" else 2e-5
    prov.RTOL = rtol
    close, at_sample = prov.close, prov.at_sample
    ok = True

    def check(name, cond, detail):
        nonlocal ok
        print(f"  {name}: {'PASS' if cond else 'FAIL'} ({detail})")
        ok &= bool(cond)

    w = a.work
    on_legs = ("prescribed_on", "serial", "low_sun", "two_stream_on", "absorbing_on", "scattering_on")
    off_legs = ("prescribed_off", "absorbing_off", "scattering_off")
    run = {leg: prov.step_dump(f"{w}/{leg}", 1) for leg in on_legs + off_legs}
    ref = run["prescribed_on"]
    same = all(set(r) == set(ref) for r in run.values())
    with_col = [leg for leg in on_legs if run[leg] and all("f_ground_lit" in f for f in run[leg].values())]
    without = [leg for leg in off_legs if run[leg] and not any("f_ground_lit" in f for f in run[leg].values())]
    walls = [k for k, f in ref.items() if k[3] != 2]
    roofs = [k for k, f in ref.items() if k[3] == 2]
    check("1. faces to compare, the new column only with the input",
          ref and same and len(with_col) == len(on_legs) and len(without) == len(off_legs) and walls and roofs,
          f"{len(ref)} faces ({len(walls)} walls, {len(roofs)} roofs), the same in {len(run)} legs: {same}; "
          f"f_ground_lit in {len(with_col)} of {len(on_legs)} legs with the input, absent in {len(without)} of "
          f"{len(off_legs)} without")
    if not ok:
        print("FAILED")
        return 1

    # The solid columns: every one has a roof face on its top.
    tops = {(k[0], k[1]): f["z_m"] for k, f in ref.items() if k[3] == 2 and k[4] < 0}
    bld = Buildings(tops, a.nx, a.ny, a.dx)
    n_ray = a.n_az * a.n_el
    count_tol = 1e-9 if a.precision == "double" else 1e-6
    summary = {}
    for leg, zen, az in (("prescribed_on", a.zenith_deg, a.azimuth_deg), ("low_sun", a.low_zenith_deg, a.low_azimuth_deg)):
        sun, cache = sun_vector(zen, az), {}
        bad, part, full, dark, ties, lit_sum, gnd_sum = [], 0, 0, 0, 0, 0.0, 0.0
        for k, f in run[leg].items():
            n_gnd, lo, hi = ground_view(f, bld, sun, a.n_az, a.n_el, cache)
            ties += hi - lo
            fg, fl = f["f_ground"], f["f_ground_lit"]
            if abs(fg - n_gnd / n_ray) > count_tol or not lo / n_ray - count_tol <= fl <= hi / n_ray + count_tol:
                bad.append((k, fg, n_gnd / n_ray, fl, lo / n_ray, hi / n_ray))
            if fg > 0.0:
                lit_sum, gnd_sum = lit_sum + fl, gnd_sum + fg
                if fl == 0.0:
                    dark += 1
                elif fl < fg:
                    part += 1
                else:
                    full += 1
        summary[leg] = (bad, part, full, dark, ties, lit_sum / gnd_sum if gnd_sum else float("nan"))
    (bad, part, full, dark, ties, share), (bad_l, part_l, full_l, dark_l, ties_l, share_l) = \
        summary["prescribed_on"], summary["low_sun"]
    check("2. the sunlit part of every face's ground view, counted from the geometry",
          not bad and not bad_l and part > 10 and full > 10 and share_l < share - 0.05 and share < 0.95
          and ties + ties_l < 0.01 * n_ray * len(ref),
          f"sun {a.zenith_deg:g} deg from the zenith: {len(bad)} of {len(ref)} faces off"
          + (f" (first {bad[0]})" if bad else "")
          + f", {full} walls see only sunlit ground, {part} both, {dark} none, sunlit share of the ground view "
          f"{share:.3f}; sun {a.low_zenith_deg:g} deg: {len(bad_l)} off"
          + (f" (first {bad_l[0]})" if bad_l else "")
          + f", {full_l} / {part_l} / {dark_l}, share {share_l:.3f}; {ties + ties_l} rays on a line between columns")

    off, on, low = run["prescribed_off"], run["prescribed_on"], run["low_sun"]
    beam_h = a.dni * math.cos(math.radians(a.zenith_deg))
    beam_l = a.dni * math.cos(math.radians(a.low_zenith_deg))
    bad_on = [k for k, f in on.items() if not close(
        f["SW_diffuse_in"], f["f_sky"] * a.diffuse + a.albedo * (f["f_ground"] * a.diffuse + f["f_ground_lit"] * beam_h))]
    bad_low = [k for k, f in low.items() if not close(
        f["SW_diffuse_in"], f["f_sky"] * a.diffuse + a.albedo * (f["f_ground"] * a.diffuse + f["f_ground_lit"] * beam_l))]
    bad_off = [k for k, f in off.items() if not close(
        f["SW_diffuse_in"], f["f_sky"] * a.diffuse + f["f_ground"] * a.albedo * (beam_h + a.diffuse))]
    bad_abs = [k for k, f in on.items() if not close(
        f["SW_abs"], (1.0 - f["albedo"]) * (f["SW_direct_in"] + f["SW_diffuse_in"]))]
    drop = max((off[k]["SW_diffuse_in"] - on[k]["SW_diffuse_in"] for k in on), default=0.0)
    check("3. the faces' own radiation: the beam reflected by the sunlit ground only, the diffuse light by all of it",
          not bad_on and not bad_low and not bad_off and not bad_abs and drop > 5.0,
          f"{len(bad_on)} faces off with the input, {len(bad_low)} at the low sun, {len(bad_off)} without the input, "
          f"{len(bad_abs)} whose absorbed light is not (1 - albedo) x incident; the input takes up to {drop:.1f} "
          f"W/m2 off a wall" + (f" (first with the input {bad_on[0]}: {on[bad_on[0]]['SW_diffuse_in']})" if bad_on else ""))

    untouched = ("SW_direct_in", "shadow", "f_sky", "f_ground", "f_bldg", "LW_in", "LW_net", "LW_ext", "T_skin", "H")
    moved = []
    for leg_on, leg_off in (("prescribed_on", "prescribed_off"), ("absorbing_on", "absorbing_off"),
                            ("scattering_on", "scattering_off")):
        for k, f in run[leg_on].items():
            g = run[leg_off][k]
            cols = untouched + (("SW_diffuse_in", "SW_abs") if k[3] == 2 else ())
            moved += [(leg_on, k, c) for c in cols if f[c] != g[c]]
    check("4. nothing else moves with the input: the beam, the shadow, the view fractions, the longwave, the roofs",
          not moved, f"{len(moved)} differences" + (f" (first {moved[0]})" if moved else ""))

    ts = run["two_stream_on"]
    s0 = a.toa * math.cos(math.radians(a.zenith_deg))
    bad_ts = [k for k, f in ts.items() if not close(f["SW_diffuse_in"], a.albedo * s0 * f["f_ground_lit"])
              or f["f_ground_lit"] != on[k]["f_ground_lit"]]
    check("5. two-stream columns, transparent sky: the reflected beam through the sunlit part of the view",
          not bad_ts and max(f["SW_diffuse_in"] for f in ts.values()) > 10.0,
          f"{len(bad_ts)} of {len(ts)} faces off" + (f" (first {bad_ts[0]}: {ts[bad_ts[0]]['SW_diffuse_in']})" if bad_ts else ""))

    ab_on, ab_off = run["absorbing_on"], run["absorbing_off"]
    bad_ab = [k for k, f in ab_on.items() if f["f_ground"] > 0.0 and not close(
        f["SW_diffuse_in"], ab_off[k]["SW_diffuse_in"] * f["f_ground_lit"] / f["f_ground"])]
    heights = sorted({k[2] for k in walls})
    check("6. absorbing sky: the light coming up is all reflected beam, shaded at every height",
          not bad_ab and len(heights) > 2 and max(f["SW_diffuse_in"] for f in ab_off.values()) > 10.0,
          f"{len(bad_ab)} of {len(walls)} walls off, in {len(heights)} cells of height"
          + (f" (first {bad_ab[0]}: {ab_on[bad_ab[0]]['SW_diffuse_in']})" if bad_ab else ""))

    sc = run["scattering_on"]
    fdir, dn, up = prov.sw_column(a.nz, a.tau, a.ssa, math.cos(math.radians(a.zenith_deg)), a.toa, a.albedo,
                                  alb_dif=a.albedo_diffuse)
    _, _, up_dark = prov.sw_column(a.nz, a.tau, a.ssa, math.cos(math.radians(a.zenith_deg)), a.toa, a.albedo,
                                   alb_dif=a.albedo_diffuse, ground_beam=0.0)
    share_beam = a.albedo * fdir[0] / up[0]
    bad_sc, gap, gap_of = [], 0.0, 0.0
    per_view = [at_sample(k, lambda m: share_beam * up[m] - (up[m] - up_dark[m])) for k in walls]
    for k, f in sc.items():
        up_f = at_sample(k, lambda m: up[m])
        expect = (f["f_sky"] * at_sample(k, lambda m: dn[m] - fdir[m])
                  + up_f * (f["f_ground"] - (f["f_ground"] - f["f_ground_lit"]) * share_beam))
        if not close(f["SW_diffuse_in"], expect):
            bad_sc.append((k, f["SW_diffuse_in"], expect))
        exact = (f["f_sky"] * at_sample(k, lambda m: dn[m] - fdir[m]) + f["f_ground_lit"] * up_f
                 + (f["f_ground"] - f["f_ground_lit"]) * at_sample(k, lambda m: up_dark[m]))
        if abs(expect - exact) > gap:
            gap, gap_of = abs(expect - exact), f["f_ground"] * up_f
    check("7. scattering sky, a ground with two albedos: the beam's share of the light coming up, shaded",
          not bad_sc and 0.3 < share_beam < 0.99 and abs(a.albedo - a.albedo_diffuse) > 0.05 and gap < a.share_gap,
          f"{len(bad_sc)} of {len(sc)} faces off the solved column"
          + (f" (first {bad_sc[0]})" if bad_sc else "")
          + f"; the beam's share of the light coming up at the ground {share_beam:.4f}; at most {gap:.3f} W/m2 "
          f"(of {gap_of:.1f} W/m2 coming up onto that face) from the column solved with a ground the beam does "
          f"not reach (allowed {a.share_gap:g}); per unit of shaded view {min(per_view):+.2f} to "
          f"{max(per_view):+.2f} W/m2 over the walls' heights, of {up[0]:.1f} W/m2 coming up at the ground")

    one_rank = run["serial"]
    differ = [k for k in on if any(on[k][c] != one_rank[k][c] for c in ("f_ground", "f_ground_lit", "SW_diffuse_in", "SW_abs"))]
    check("8. one rank and several: the same sunlit parts and the same light", not differ,
          f"{len(differ)} of {len(on)} faces differ" + (f" (first {differ[0]})" if differ else ""))

    cont, rest = prov.step_dump(f"{w}/prescribed_on", 2), prov.step_dump(f"{w}/restart", 2)
    differ = [k for k in cont if k not in rest or any(cont[k][c] != rest[k][c] for c in cont[k])]
    check("9. a restart finds the sunlit parts again", cont and set(cont) == set(rest) and not differ
          and any(f["f_ground_lit"] > 0.0 for f in rest.values()),
          f"{len(cont)} faces at step 2, {len(differ)} differ after the restart"
          + (f" (first {differ[0]})" if differ else ""))

    good, detail = True, []
    for lev, name in ((0, "faces/set"), (1, "faces/set.lev1")):
        t_on, t_off = prov.faces(f"{w}/two_level_on/{name}.rank*.csv"), prov.faces(f"{w}/two_level_off/{name}.rank*.csv")
        if not t_on or set(t_on) != set(t_off) or not all("f_ground_lit" in f for f in t_on.values()):
            good = False
            detail.append(f"level {lev}: {len(t_on)} and {len(t_off)} faces, or no f_ground_lit")
            continue
        over = [k for k, f in t_on.items() if f["f_ground_lit"] > f["f_ground"] or f["f_ground_lit"] < 0.0]
        bad = [k for k, f in t_on.items() if f["f_ground"] > 0.0 and not close(
            f["SW_diffuse_in"], t_off[k]["SW_diffuse_in"] * f["f_ground_lit"] / f["f_ground"])]
        roofs_moved = [k for k, f in t_on.items() if k[3] == 2 and f["SW_abs"] != t_off[k]["SW_abs"]]
        seeing = [f for f in t_on.values() if f["f_ground"] > 0.0]
        part = [f for f in seeing if 0.0 < f["f_ground_lit"] < f["f_ground"]]
        dark = [f for f in seeing if f["f_ground_lit"] == 0.0]
        lit_share = sum(f["f_ground_lit"] for f in seeing) / sum(f["f_ground"] for f in seeing) if seeing else float("nan")
        good &= not over and not bad and not roofs_moved and len(part) > 5 and lit_share < 0.95
        detail.append(f"level {lev}: {len(t_on)} faces, {len(seeing)} seeing the ground ({len(part)} partly shaded, "
                      f"{len(dark)} all shaded), sunlit share {lit_share:.3f}, {len(over)} out of range, {len(bad)} off, "
                      f"{len(roofs_moved)} roofs moved")
    check("10. two levels: each level's faces, their sunlit parts and the reflected light", good, "; ".join(detail))

    print("ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
