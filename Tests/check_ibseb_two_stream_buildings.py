#!/usr/bin/env python3
"""Check the two-stream columns and ground balance against the buildings.

    python3 check_ibseb_two_stream_buildings.py <work dir> --cosz C --zenith-deg Z --azimuth-deg A
            --toa S0 --albedo a --tau t --nz N [--two-level-cosz C1 --two-level-zenith-deg Z1
            --two-level-azimuth-deg A1] [--precision double|single]

Reads the legs Tests/RunIBSEBTwoStreamBuildings.cmake ran (one directory each) with
erf.ibseb.two_stream_buildings: the 2D plotfiles (seb_sw_net, ibseb_ground_shadow,
seb_t_sfc) and the 3D plotfiles (terrain_IB_mask, qsrc_sw, qsrc_lw). The buildings are
taken from the blanking of the 3D plotfile (terrain_IB_mask >= 0.5); the shadow they cast
is computed here on its own, by intersecting the ray from the centre of each column's
ground with every solid column as a box (slab method, the periodic images included),
not by the code's column walk. Each check names the defect it guards; the exit code is
the number of failed checks.

  1. clear: the ground shadow is that geometry's, column for column, on a sun from the
     south-south-west (a wrong sun direction or a ray from the wrong point fails; the
     periodic wrap is covered by the unit test IBSEBGroundShadow);
  2. clear: the net shortwave of the ground balance is (1 - albedo) S0 cos z on the sunlit
     open ground, zero on the shaded open ground (the sky has no diffuse light) and zero
     under the buildings (a footprint or a shadow left in the sun fails);
  3. clear: under a building the ground keeps its initial skin temperature after two
     steps, and the sunlit open ground does not (a footprint the force-restore update
     still advances fails); the mean skin temperature the radiation diagnostics report
     after each update is the mean over the open ground (one that counts the footprints
     fails);
  4. absorbing: the sunlit open ground takes (1 - albedo) S0 cos z exp(-tau N / cos z), the
     beam the sweep brings through N absorbing layers, and the shaded ground none (a beam
     taken from the wrong interface, or the total shortwave taken off, fails);
  5. absorbing: no two-stream heating in a solid cell, while the air at the same heights
     is heated (a mask left off, or applied to the air, fails);
  6. scattering: the sunlit and the shaded open ground differ by exactly the absorbed
     beam, (1 - albedo) S0 cos z exp(-tau N / cos z), so the shaded ground keeps the sky's
     diffuse light (taking off the total shortwave fails), while the ground under the
     buildings gets none of it (a footprint given the shaded ground's light fails: the
     footprints are in shadow too, so only a sky with diffuse light tells them apart);
  7. off (erf.ibseb.two_stream_buildings = false): the columns are blind to the buildings
     as before the input: the sunlit value everywhere, the footprints included, no shadow
     field, and heating in the solid cells;
  8. restart: two steps straight and one step plus a restart give the same ground fields
     bit for bit (a shadow or a held footprint lost on restart fails); the plotfile a run
     writes before its first step, and the one a restart writes, carry the shadow of their
     sun (the step-0 one the geometry's, the restart's that of the straight run at that
     step) and -999 for the net shortwave, which no update has made yet (a field left at
     its start, or a shadow of zero, fails);
  9. night: the sun below the horizon shades every column and the ground gets no
     shortwave (with the traps on);
 10. two levels: on both levels the shadow is the geometry's, with the tower outside the
     refined level shading level 1 through the coarse column map (a refined level that
     ignores the buildings outside it fails);
 11. a refined patch below the domain top, whose heating rates come from level 0 by
     interpolation: no heating in its solid cells, heating in its air at those heights
     whose level-0 cell is air too (the interpolation fills the building's cells, and a
     patch left unmasked fails; an air cell of the patch inside a solid cell of level 0,
     where level 0's building fills the cell beyond its footprint, takes its heating from
     level 0's zero by the limited interpolation and is left out);
 12. fixed shortwave (erf.radiation.seb_use_radiation_fluxes = false): the footprints are
     held and get no shortwave, the open ground gets erf.radiation.seb_sw_flux_default,
     no shadow is cast (the plot variable is -999), and ERF says so at start-up.
"""
import argparse
import math
import os
import re
import struct
import sys


# ---------------------------------------------------------------------------
# Plotfiles (standard library only): a field is a dict {(i, j[, k]): value} over the
# cells a level has.
# ---------------------------------------------------------------------------

def read_plotfile(path):
    """(names, [boxes per level], dx per level) with each box (lo, (ni, nj, nk), {name: flat values})."""
    with open(os.path.join(path, "Header")) as f:
        lines = f.read().split("\n")
    nvar = int(lines[1])
    names = lines[2:2 + nvar]
    p = 2 + nvar + 2
    finest = int(lines[p]); p += 1
    p += 3 + 2
    dxs = []
    for _ in range(finest + 1):
        dxs.append([float(v) for v in lines[p].split()]); p += 1
    levels = []
    for lev in range(finest + 1):
        boxes = []
        with open(os.path.join(path, f"Level_{lev}", "Cell_H")) as f:
            ch = f.read()
        for name, off in re.findall(r"FabOnDisk: (\S+) (\d+)", ch):
            with open(os.path.join(path, f"Level_{lev}", name), "rb") as f:
                f.seek(int(off))
                head = b""
                while not head.endswith(b"\n"):
                    head += f.read(1)
                m = re.search(rb"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\) \(\d+,\d+,\d+\)\) (\d+)", head)
                lo = [int(m.group(k)) for k in (1, 2, 3)]
                hi = [int(m.group(k)) for k in (4, 5, 6)]
                nc = int(m.group(7))
                n = [hi[d] - lo[d] + 1 for d in range(3)]
                little = b"(8 7 6 5 4 3 2 1)" in head or b"(4 3 2 1)" in head
                size = 8 if b"(8 7 6 5 4 3 2 1)" in head or b"(1 2 3 4 5 6 7 8)" in head else 4
                count = n[0] * n[1] * n[2]
                fmt = ("<" if little else ">") + ("d" if size == 8 else "f") * count
                data = {}
                for c in range(nc):
                    data[names[c]] = struct.unpack(fmt, f.read(size * count))
                boxes.append((lo, n, data))
        levels.append(boxes)
    return names, levels, dxs


def field(levels, lev, name, two_d=False):
    """{(i, j): v} for a 2D plotfile (its k = 0 slab), {(i, j, k): v} otherwise."""
    out = {}
    for lo, n, data in levels[lev]:
        vals = data[name]
        for k in range(n[2]):
            for j in range(n[1]):
                for i in range(n[0]):
                    v = vals[(k * n[1] + j) * n[0] + i]
                    if two_d:
                        out[(lo[0] + i, lo[1] + j)] = v
                    else:
                        out[(lo[0] + i, lo[1] + j, lo[2] + k)] = v
    return out


def plotfiles(run, prefix):
    pat = re.compile(r"^" + re.escape(prefix) + r"(\d+)$")
    found = []
    for name in os.listdir(run):
        m = pat.match(name)
        if m:
            found.append((int(m.group(1)), os.path.join(run, name)))
    if not found:
        raise FileNotFoundError(f"no {prefix}* plotfile in {run}")
    return sorted(found)


def plotfile(run, prefix, step=None):
    found = plotfiles(run, prefix)
    if step is None:
        return found[-1][1]
    for s_, p in found:
        if s_ == step:
            return p
    raise FileNotFoundError(f"no {prefix}{step:05d} in {run}")


def plotfile_time(path):
    with open(os.path.join(path, "Header")) as f:
        lines = f.read().split("\n")
    return float(lines[2 + int(lines[1]) + 1])


# ---------------------------------------------------------------------------
# The buildings' shadow, on its own
# ---------------------------------------------------------------------------

def sun_vector(zenith_deg, azimuth_deg):
    """Unit vector toward the sun: azimuth clockwise from north, x east, y north."""
    z, a = math.radians(zenith_deg), math.radians(azimuth_deg)
    return math.sin(a) * math.sin(z), math.cos(a) * math.sin(z), math.cos(z)


def tops_from_mask(mask, dz):
    """{(i, j): top height} of the solid columns (terrain_IB_mask >= 0.5)."""
    tops = {}
    for (i, j, k), v in mask.items():
        if v >= 0.5:
            tops[(i, j)] = max(tops.get((i, j), 0.0), (k + 1) * dz)
    return tops


def geometric_shadow(columns, tops, nx, ny, dx, dy, s, periodic=(True, True)):
    """{column: (shadow, grazing)} for the given (i, j) columns: the ray from the centre of
    the column's ground along s meets the box [x, x + dx] x [y, y + dy] x [0, top] of a
    solid column, or of a periodic image, below its top on entry (slab method); grazing
    when it passes within a grazing band of a box's top or crosses one in less than that band
    (1e-4 of the column width, above the code's 1e-6 band and single precision's resolution)."""
    sx, sy, sz = s
    band = 1.0e-4 * min(dx, dy)
    out = {}
    shx = [-nx * dx, 0.0, nx * dx] if periodic[0] else [0.0]
    shy = [-ny * dy, 0.0, ny * dy] if periodic[1] else [0.0]
    big = 1.0e300
    for (i, j) in columns:
        if sz <= 0.0 or (i, j) in tops:
            out[(i, j)] = (1, False)
            continue
        xc, yc = (i + 0.5) * dx, (j + 0.5) * dy
        hit = 0
        graze = False
        for (ci, cj), top in tops.items():
            for ox in shx:
                for oy in shy:
                    x0, x1 = ci * dx + ox, (ci + 1) * dx + ox
                    y0, y1 = cj * dy + oy, (cj + 1) * dy + oy
                    if sx != 0.0:
                        ta, tb = (x0 - xc) / sx, (x1 - xc) / sx
                        tx0, tx1 = min(ta, tb), max(ta, tb)
                    elif x0 < xc < x1:
                        tx0, tx1 = -big, big
                    else:
                        continue
                    if sy != 0.0:
                        ta, tb = (y0 - yc) / sy, (y1 - yc) / sy
                        ty0, ty1 = min(ta, tb), max(ta, tb)
                    elif y0 < yc < y1:
                        ty0, ty1 = -big, big
                    else:
                        continue
                    t_in = max(tx0, ty0, 0.0)
                    t_out = min(tx1, ty1)
                    if t_out > t_in - band and (abs(t_in * sz - top) < band or abs(t_out - t_in) < band):
                        graze = True
                    if t_out > t_in and t_in * sz < top:
                        hit = 1
        out[(i, j)] = (hit, graze)
    return out


class Checks:
    def __init__(self):
        self.failed = 0

    def __call__(self, name, ok, detail):
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")
        if not ok:
            self.failed += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("--cosz", type=float, required=True)
    ap.add_argument("--zenith-deg", type=float, required=True)
    ap.add_argument("--azimuth-deg", type=float, required=True)
    ap.add_argument("--toa", type=float, required=True)
    ap.add_argument("--albedo", type=float, required=True)
    ap.add_argument("--tau", type=float, required=True)
    ap.add_argument("--nz", type=int, required=True)
    ap.add_argument("--two-level-zenith-deg", type=float, required=True)
    ap.add_argument("--two-level-azimuth-deg", type=float, required=True)
    ap.add_argument("--precision", default="double")
    a = ap.parse_args()
    rtol = 1e-10 if a.precision == "double" else 2e-5
    check = Checks()
    s = sun_vector(a.zenith_deg, a.azimuth_deg)
    if abs(s[2] - a.cosz) > 1e-12:
        print(f"FAIL  the sun: cos z of zenith {a.zenith_deg} deg is {s[2]}, not {a.cosz}")
        return 1
    clear_sun = (1.0 - a.albedo) * a.toa * a.cosz
    beam_abs = (1.0 - a.albedo) * a.toa * a.cosz * math.exp(-a.tau * a.nz / a.cosz)

    def leg(name):
        return os.path.join(a.work, name)

    def ground(run, name, step=None, lev=0):
        names, levels, _ = read_plotfile(plotfile(run, "plt2d", step))
        return field(levels, lev, name, two_d=True)

    def volume(run, name, step, lev=0):
        names, levels, dxs = read_plotfile(plotfile(run, "plt", step))
        return field(levels, lev, name), dxs[lev]

    def rel(x, ref):
        return abs(x - ref) / abs(ref)

    # ---- clear -------------------------------------------------------------
    mask, dx = volume(leg("clear"), "terrain_IB_mask", 0)
    nx = max(i for i, _, _ in mask) + 1
    ny = max(j for _, j, _ in mask) + 1
    tops = tops_from_mask(mask, dx[2])
    foot = {c for c in ((i, j) for i in range(nx) for j in range(ny)) if mask[(c[0], c[1], 0)] >= 0.5}
    columns = [(i, j) for i in range(nx) for j in range(ny)]
    geo = geometric_shadow(columns, tops, nx, ny, dx[0], dx[1], s)
    shadow = ground(leg("clear"), "ibseb_ground_shadow")
    sw = ground(leg("clear"), "seb_sw_net")
    near = {c for c in columns if geo[c][1]}
    lit = [c for c in columns if c not in foot and c not in near and geo[c][0] == 0]
    shaded = [c for c in columns if c not in foot and c not in near and geo[c][0] == 1]
    mism = [c for c in columns if c not in near and shadow[c] != geo[c][0]]
    check("1. clear: the ground shadow is the buildings' geometry",
          not mism and shaded and lit,
          f"{len(shaded)} shaded and {len(lit)} sunlit open columns, {len(foot)} under the buildings, "
          f"{len(near)} grazing left out; {len(mism)} differ from the geometry")
    err_lit = max(rel(sw[c], clear_sun) for c in lit)
    check("2. clear: net shortwave: sunlit (1 - a) S0 cos z, shaded and footprints zero",
          err_lit <= rtol and all(sw[c] == 0.0 for c in shaded) and all(sw[c] == 0.0 for c in foot),
          f"sunlit {min(sw[c] for c in lit):.6f} to {max(sw[c] for c in lit):.6f} W/m2 (expected {clear_sun:.6f}, "
          f"rel. error {err_lit:.1e}); shaded max {max(abs(sw[c]) for c in shaded):.3g}; "
          f"footprints max {max(abs(sw[c]) for c in foot):.3g}")
    t0 = ground(leg("clear"), "seb_t_sfc", 0)
    t2 = ground(leg("clear"), "seb_t_sfc")
    held = all(t2[c] == t0[c] for c in foot)
    moved = min(abs(t2[c] - t0[c]) for c in lit)
    # The diagnostics' mean skin after each update, against the plotfile's open ground.
    with open(os.path.join(leg("clear"), "radiation_diag.csv")) as f:
        head = f.readline().strip().split(",")
        rows = [dict(zip(head, line.strip().split(","))) for line in f if line.strip()]
    diag_err, n_rows = 0.0, 0
    open_cols = [c for c in columns if c not in foot]
    for st in (1, 2):
        t_st = ground(leg("clear"), "seb_t_sfc", st)
        t_plot = plotfile_time(plotfile(leg("clear"), "plt2d", st))
        post = [r for r in rows if r["call_site"] == "post_dycore" and r["level"] == "0"
                and abs(float(r["time"]) - t_plot) < 1e-9]
        if post:
            n_rows += 1
            mean_open = sum(t_st[c] for c in open_cols) / len(open_cols)
            diag_err = max(diag_err, rel(float(post[-1]["T_s_mean"]), mean_open))
    counted = len(foot) / len(columns) * moved / 300.0
    check("3. clear: the ground under the buildings is held, the sunlit ground is not; the means leave the footprints out",
          held and moved > 0.0 and n_rows == 2 and diag_err <= 5e-6,
          f"footprints {max(abs(t2[c] - t0[c]) for c in foot):.3g} K from their start; sunlit ground at least {moved:.3g} K; "
          f"reported mean skin vs the open ground's in {n_rows} updates: rel. error {diag_err:.1e} "
          f"(the CSV holds 7 digits; counting the footprints would give about {counted:.1e})")

    # ---- absorbing -----------------------------------------------------------
    sw = ground(leg("absorbing"), "seb_sw_net")
    err = max(rel(sw[c], beam_abs) for c in lit)
    check("4. absorbing: the sunlit ground takes the beam through the column, the shaded none",
          err <= rtol and all(sw[c] == 0.0 for c in shaded) and all(sw[c] == 0.0 for c in foot),
          f"sunlit {min(sw[c] for c in lit):.6f} to {max(sw[c] for c in lit):.6f} W/m2 (expected {beam_abs:.6f}, "
          f"rel. error {err:.1e}); shaded max {max(abs(sw[c]) for c in shaded):.3g}")
    qsw, _ = volume(leg("absorbing"), "qsrc_sw", 1)
    qlw, _ = volume(leg("absorbing"), "qsrc_lw", 1)
    solid = [c for c, v in mask.items() if v >= 0.5]
    k_bld = {c[2] for c in solid}
    air = [c for c, v in mask.items() if v < 0.5 and c[2] in k_bld]
    check("5. absorbing: no heating in the solid cells, heating in the air beside them",
          all(qsw[c] == 0.0 and qlw[c] == 0.0 for c in solid)
          and min(abs(qsw[c]) for c in air) > 0.0 and min(abs(qlw[c]) for c in air) > 0.0,
          f"{len(solid)} solid cells: |SW| max {max(abs(qsw[c]) for c in solid):.3g}, "
          f"|LW| max {max(abs(qlw[c]) for c in solid):.3g} K/s; air at those heights: "
          f"|SW| {min(abs(qsw[c]) for c in air):.3g} to {max(abs(qsw[c]) for c in air):.3g}, "
          f"|LW| {min(abs(qlw[c]) for c in air):.3g} to {max(abs(qlw[c]) for c in air):.3g} K/s")

    # ---- scattering ----------------------------------------------------------
    sw = ground(leg("scattering"), "seb_sw_net")
    v_lit = [sw[c] for c in lit]
    v_sh = [sw[c] for c in shaded]
    d_lit, d_sh = sum(v_lit) / len(v_lit), sum(v_sh) / len(v_sh)
    spread = max(max(v_lit) - min(v_lit), max(v_sh) - min(v_sh))
    err = rel(d_lit - d_sh, beam_abs)
    foot_max = max(abs(sw[c]) for c in foot)
    check("6. scattering: sunlit minus shaded is the absorbed beam; the shaded ground keeps the diffuse light",
          err <= 1e3 * rtol and d_sh > 0.0 and spread <= 1e3 * rtol * d_lit and foot_max == 0.0,
          f"sunlit {d_lit:.6f}, shaded {d_sh:.6f} W/m2 (spread {spread:.2g}); difference {d_lit - d_sh:.6f}, "
          f"expected {beam_abs:.6f} (rel. error {err:.1e}); footprints max {foot_max:.3g}")

    # ---- off -----------------------------------------------------------------
    sw = ground(leg("off"), "seb_sw_net")
    sh = ground(leg("off"), "ibseb_ground_shadow")
    qsw_off, _ = volume(leg("off"), "qsrc_sw", 1)
    err = max(rel(v, beam_abs) for v in sw.values())
    check("7. off: the columns are blind to the buildings, as before the input",
          err <= rtol and all(v == -999.0 for v in sh.values()) and min(abs(qsw_off[c]) for c in solid) > 0.0,
          f"net shortwave {min(sw.values()):.6f} to {max(sw.values()):.6f} W/m2 everywhere, footprints included "
          f"(rel. error {err:.1e}); shadow field {sorted(set(sh.values()))}; solid cells heated, "
          f"|SW| min {min(abs(qsw_off[c]) for c in solid):.3g} K/s")

    # ---- restart -------------------------------------------------------------
    same, detail, n_compared = True, [], 0
    for prefix in ("plt2d", "plt"):
        _, la, _ = read_plotfile(plotfile(leg("clear"), prefix, 2))
        _, lb, _ = read_plotfile(plotfile(leg("restart"), prefix, 2))
        for v in ("seb_sw_net", "ibseb_ground_shadow", "seb_t_sfc", "qsrc_sw", "qsrc_lw", "theta"):
            if v not in la[0][0][2]:
                continue
            ok = field(la, 0, v) == field(lb, 0, v)
            same &= ok
            n_compared += 1
            detail.append(f"{v} {'=' if ok else 'differs'}")
    # The plotfiles written before any update: the step-0 one and the restart's.
    sh0 = ground(leg("clear"), "ibseb_ground_shadow", 0)
    sw0 = ground(leg("clear"), "seb_sw_net", 0)
    sh_r = ground(leg("restart"), "ibseb_ground_shadow", 1)
    sw_r = ground(leg("restart"), "seb_sw_net", 1)
    sh_1 = ground(leg("clear"), "ibseb_ground_shadow", 1)
    first_ok = (all(sh0[c] == geo[c][0] for c in columns if c not in near) and sh_r == sh_1
                and all(v == -999.0 for v in sw0.values()) and all(v == -999.0 for v in sw_r.values()))
    check("8. restart: one step and a restart give the two-step ground fields bit for bit; the first plotfiles",
          same and n_compared == 6 and first_ok,
          ", ".join(detail) + f" ({n_compared} of 6 compared); step 0: {int(sum(sh0.values()))} shaded columns, "
          f"net shortwave {sorted(set(sw0.values()))}; the restart's plotfile: shadow "
          f"{'as the straight run' if sh_r == sh_1 else 'differs from the straight run'}, net shortwave {sorted(set(sw_r.values()))}")

    # ---- night ---------------------------------------------------------------
    sw = ground(leg("night"), "seb_sw_net")
    sh = ground(leg("night"), "ibseb_ground_shadow")
    check("9. night: every column shaded and no shortwave on the ground",
          all(v == 1.0 for v in sh.values()) and all(v == 0.0 for v in sw.values()),
          f"shadow {sorted(set(sh.values()))}, net shortwave {min(sw.values()):.3g} to {max(sw.values()):.3g} W/m2")

    # ---- two levels ----------------------------------------------------------
    s2 = sun_vector(a.two_level_zenith_deg, a.two_level_azimuth_deg)
    m0, dx0 = volume(leg("two_level"), "terrain_IB_mask", 0, 0)
    m1, dx1 = volume(leg("two_level"), "terrain_IB_mask", 0, 1)
    nx0 = max(i for i, _, _ in m0) + 1
    ny0 = max(j for _, j, _ in m0) + 1
    tops0 = tops_from_mask(m0, dx0[2])
    tops1_own = tops_from_mask(m1, dx1[2])
    on1 = {(i, j) for (i, j, k) in m1}
    # Level 1 sees its own buildings where it has cells and level 0's elsewhere.
    tops1 = dict(tops1_own)
    for (i, j), top in tops0.items():
        for di in (0, 1):
            for dj in (0, 1):
                c = (2 * i + di, 2 * j + dj)
                if c not in on1:
                    tops1[c] = top
    out, ok_all = [], True
    for lev, tops_l, dxl, nxl, nyl, cols in ((0, tops0, dx0, nx0, ny0, [(i, j) for i in range(nx0) for j in range(ny0)]),
                                             (1, tops1, dx1, 2 * nx0, 2 * ny0, sorted(on1))):
        geo_l = geometric_shadow(cols, tops_l, nxl, nyl, dxl[0], dxl[1], s2)
        sh = ground(leg("two_level"), "ibseb_ground_shadow", None, lev)
        have = [c for c in cols if not geo_l[c][1]]
        bad = [c for c in have if sh[c] != geo_l[c][0]]
        n_sh = sum(geo_l[c][0] for c in have)
        ok_all &= (not bad) and n_sh > 0
        out.append(f"level {lev}: {len(have)} columns, {n_sh} shaded, {len(bad)} differ")
    own = geometric_shadow(sorted(on1), tops1_own, 2 * nx0, 2 * ny0, dx1[0], dx1[1], s2)
    full = geometric_shadow(sorted(on1), tops1, 2 * nx0, 2 * ny0, dx1[0], dx1[1], s2)
    from_outside = sum(1 for c in on1 if full[c][0] == 1 and own[c][0] == 0)
    check("10. two levels: the shadow is the geometry's on both levels, the outside tower's on level 1 too",
          ok_all and from_outside > 0,
          "; ".join(out) + f"; {from_outside} columns of level 1 shaded only by buildings outside it")

    # ---- a nested patch --------------------------------------------------------
    m0n, dxn0 = volume(leg("nested"), "terrain_IB_mask", 1, 0)
    m1, dxn1 = volume(leg("nested"), "terrain_IB_mask", 1, 1)
    rr = [int(round(dxn0[d] / dxn1[d])) for d in range(3)]
    q1sw, _ = volume(leg("nested"), "qsrc_sw", 1, 1)
    q1lw, _ = volume(leg("nested"), "qsrc_lw", 1, 1)
    solid1 = [c for c, v in m1.items() if v >= 0.5]
    k1 = {c[2] for c in solid1}
    air1 = [c for c, v in m1.items() if v < 0.5 and c[2] in k1 and m0n[(c[0] // rr[0], c[1] // rr[1], c[2] // rr[2])] < 0.5]
    check("11. a patch below the top, its rates interpolated: no heating in its solid cells, heating in its air",
          solid1 and all(q1sw[c] == 0.0 and q1lw[c] == 0.0 for c in solid1)
          and min(abs(q1sw[c]) for c in air1) > 0.0 and min(abs(q1lw[c]) for c in air1) > 0.0,
          f"level 1 up to k = {max(c[2] for c in m1)}: {len(solid1)} solid cells, |SW| max "
          f"{max(abs(q1sw[c]) for c in solid1):.3g}, |LW| max {max(abs(q1lw[c]) for c in solid1):.3g} K/s; "
          f"air at those heights |SW| {min(abs(q1sw[c]) for c in air1):.3g} to {max(abs(q1sw[c]) for c in air1):.3g} K/s")

    # ---- fixed shortwave -------------------------------------------------------
    sw = ground(leg("fixed_sw"), "seb_sw_net")
    sh = ground(leg("fixed_sw"), "ibseb_ground_shadow")
    tf0 = ground(leg("fixed_sw"), "seb_t_sfc", 0)
    tf2 = ground(leg("fixed_sw"), "seb_t_sfc")
    with open(os.path.join(leg("fixed_sw"), "simulation.log")) as f:
        noted = "takes the fixed erf.radiation.seb_sw_flux_default" in f.read()
    check("12. fixed shortwave: footprints held with no shortwave, the open ground the default, no shadow cast",
          noted and all(v == -999.0 for v in sh.values()) and all(sw[c] == 0.0 for c in foot)
          and all(sw[c] == 300.0 for c in open_cols) and all(tf2[c] == tf0[c] for c in foot),
          f"note printed: {noted}; shadow field {sorted(set(sh.values()))}; net shortwave under the buildings "
          f"{sorted(set(sw[c] for c in foot))}, open {sorted(set(sw[c] for c in open_cols))}; footprints "
          f"{max(abs(tf2[c] - tf0[c]) for c in foot):.3g} K from their start")

    return check.failed


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, KeyError, ValueError, IndexError) as err:
        # A leg that did not run, or a plotfile without a field the checks read.
        print(f"FAIL  the runs could not be read: {type(err).__name__}: {err}")
        sys.exit(1)
