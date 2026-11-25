"""Built-in tests (``erfdomainadvisor.py --self-test``), all offline.

Each test checks a number against a closed form or a published value, not against the code's
own output: the HRRR grid corners, the level chain, the geometric vertical grid, the slope of a
Gaussian hill, a PNG round trip through every filter, and the box arithmetic ERF uses to place a
nest from its met_em file.  Each test was seen to fail against a one-line mutant of the code it
covers (see README.md).
"""

import math
import os
import random
import struct
import subprocess
import sys
import tempfile
import zlib

from . import writers
from .phenomena import PHENOMENA
from .plan import (ERF_FIXED_DT_CFL, NICE_DT, SOUND_SPEED, PlanError, Settings, base_time_step,
                   build, erf_nest_box, fire_grid, level_chain, regime, validate, vertical_grid)
from .projection import (HRRR, HRRR_DX, HRRR_NX, HRRR_NY, Lambert, hrrr_index)
from .terrain import (AsciiGridDEM, AnalyticDEM, TerrainError, TerrariumDEM, boundary_mismatch,
                      cell_slope, check, decode_png, despike, slope_stats, level_terrain)

TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _close(a, b, tol, what):
    assert abs(a - b) <= tol, "%s: %r, expected %r (tolerance %g)" % (what, a, b, tol)


def _raises(func, text):
    try:
        func()
    except PlanError as exc:
        assert text in str(exc), "message %r does not name %r" % (str(exc), text)
        return
    raise AssertionError("no PlanError naming %r" % text)


# ---- projection ---------------------------------------------------------------------------

def test_hrrr_corners():
    # First point and far corner of the HRRR CONUS grid (NCEP grid 184 metadata).
    i, j = hrrr_index(21.138123, -122.719528)
    _close(i, 0.0, 1e-3, "i of the first point")
    _close(j, 0.0, 1e-3, "j of the first point")
    i, j = hrrr_index(47.842195, -60.917193)
    _close(i, HRRR_NX - 1, 0.02, "i of the far corner")
    _close(j, HRRR_NY - 1, 0.02, "j of the far corner")


def test_lambert_round_trip():
    rnd = random.Random(1)
    proj = Lambert(38.5, -97.5)
    for _ in range(50):
        lat, lon = rnd.uniform(20, 55), rnd.uniform(-130, -60)
        la, lo = proj.inverse(*proj.forward(lat, lon))
        _close(la, lat, 1e-9, "latitude round trip")
        _close(lo, lon, 1e-9, "longitude round trip")
    # map factor is 1 on the standard parallel
    _close(HRRR.map_factor(38.5), 1.0, 1e-12, "map factor at truelat")


# ---- levels, vertical grid, regimes -------------------------------------------------------

def test_level_chain():
    n, b = level_chain(900.0 / 27, 3, 1000.0)
    assert n == 3, n
    _close(b, 900.0, 1e-9, "base for 33.3 m, ratio 3")
    n, b = level_chain(100.0, 2, 1000.0)
    assert n == 3, n
    _close(b, 800.0, 1e-9, "base for 100 m, ratio 2")
    assert level_chain(2000.0, 3, 1000.0) == (0, 2000.0)
    _raises(lambda: level_chain(1.0, 2, 3000.0), "at most 6")


def test_validation_aborts():
    base = dict(lat=40.0, lon=-105.0, phenomenon="wildfire")
    for kw, text in ((dict(lat=95.0), "--lat"), (dict(lat=-90.0), "--lat"), (dict(ratio=4), "--ratio"),
                     (dict(target_dx=0.0), "--target-dx"), (dict(phenomenon="tornado"), "--phenomenon"),
                     (dict(fire_dx=0.1), "--fire-dx"), (dict(max_stretch=1.0), "--max-stretch"),
                     (dict(nest_margin=2), "--nest-margin")):
        args = dict(base)
        args.update(kw)
        _raises(lambda: validate(Settings(**args)), text)


def test_vertical_grid():
    for top, dz0, dz_top, qmax in ((10000.0, 10.0, 500.0, 1.08), (20000.0, 50.0, 600.0, 1.08),
                                   (4000.0, 2.0, 400.0, 1.05), (6000.0, 20.0, 400.0, 1.08)):
        nz, q, dzt = vertical_grid(top, dz0, dz_top, qmax)
        assert nz % 8 == 0, "nz %d is not a multiple of 8" % nz
        _close(dz0 * (q ** nz - 1) / (q - 1), top, 1e-6 * top, "sum of the cell heights")
        assert q <= qmax and dzt <= dz_top + 1e-9, (q, dzt)
        # the next smaller multiple of 8 must break a rule
        if nz > 8:
            from .plan import solve_stretch
            q8 = solve_stretch(nz - 8, dz0, top)
            assert q8 > qmax or dz0 * q8 ** (nz - 9) > dz_top, "nz %d is not the smallest" % nz


def test_regime():
    assert regime(1500.0, 1500.0) == "PBL scheme"
    assert regime(299.0, 1500.0) == "LES"
    assert regime(300.0, 1500.0) == "gray zone"     # dx = 0.2 zi is not LES
    assert regime(600.0, 1500.0) == "gray zone"


# ---- nests --------------------------------------------------------------------------------

def _plan(**kw):
    args = dict(lat=34.07, lon=-118.55, phenomenon="wildfire", wind_from=45.0)
    args.update(kw)
    return build(Settings(**args))


def test_nests_and_erf_box():
    for kw in ({}, dict(phenomenon="convective_bl_les", wind_from=None, lat=36.6, lon=-97.5),
               dict(phenomenon="downslope_windstorm", ratio=2), dict(phenomenon="urban_canopy")):
        p = _plan(**kw)
        r = p.ratio
        for lv in p.levels[1:]:
            par = p.levels[lv.lev - 1]
            assert lv.nx % r == 0 and lv.ny % r == 0, "level %d size not a multiple of %d" % (lv.lev, r)
            assert lv.i_start + lv.nx // r <= par.nx and lv.j_start + lv.ny // r <= par.ny
            need = p.settings.nest_margin if lv.lev >= 2 else 1
            assert lv.margin >= need, "level %d margin %d < %d" % (lv.lev, lv.margin, need)
            _close(lv.dx * r, par.dx, 1e-9, "dx ratio")
            _close(lv.x1 - lv.x0, lv.nx * lv.dx, 1e-6, "level %d x extent vs nx dx" % lv.lev)
            _close(lv.y1 - lv.y0, lv.ny * lv.dx, 1e-6, "level %d y extent vs ny dx" % lv.lev)
            # ERF's box from the met_em attributes holds the point of interest
            ilo, jlo, ihi, jhi = erf_nest_box(p.levels, r, lv.lev)
            ip, jp = p.point_xy[0] / lv.dx, p.point_xy[1] / lv.dx
            assert ilo <= ip <= ihi + 1 and jlo <= jp <= jhi + 1, \
                "point (%g, %g) outside ERF box %s at level %d" % (ip, jp, (ilo, jlo, ihi, jhi), lv.lev)
            _close(ilo * lv.dx, lv.x0, 1e-6, "ERF box x0 vs plan")
            _close(jlo * lv.dx, lv.y0, 1e-6, "ERF box y0 vs plan")


def test_namelist_matches_plan():
    p = _plan(start="2025-01-07_12")
    nl = writers.namelist_wps(p)

    def row(key):
        line = [ln for ln in nl if ln.strip().startswith(key + " ")][0]
        return [int(v) for v in line.split("=")[1].replace(",", " ").split()]
    assert row("e_we") == [lv.nx + 1 for lv in p.levels]
    assert row("e_sn") == [lv.ny + 1 for lv in p.levels]
    assert row("i_parent_start")[1:] == [lv.i_start + 1 for lv in p.levels[1:]]
    assert row("j_parent_start")[1:] == [lv.j_start + 1 for lv in p.levels[1:]]
    assert all((e - 1) % p.ratio == 0 for e in row("e_we")[1:])


def test_fire_grid():
    assert fire_grid(900.0 / 27, 10.0, 3) == (4, 60), fire_grid(900.0 / 27, 10.0, 3)
    assert fire_grid(100.0, 10.0, 3)[0] == 10
    assert fire_grid(5.0, 10.0, 3)[0] == 1
    _raises(lambda: fire_grid(1000.0, 10.0, 3), "--fire-dx")
    p = _plan()
    assert p.fire and PHENOMENA["wildfire"].fire
    fl = p.levels[-1]
    assert p.fire_box % p.fire_ratio == 0 and p.fire_box % p.ratio == 0
    assert fl.nx % p.fire_box == 0 and fl.ny % p.fire_box == 0, (fl.nx, fl.ny, p.fire_box)
    assert fl.margin >= p.settings.nest_margin, "fire level margin %d" % fl.margin
    deck = writers.inputs_deck(p, "test")
    assert "erf.fire.grid_ratio   = 4      # fire cells of 8.33333 m (--fire-dx 10 m)" in deck
    mgs = [ln for ln in deck if ln.startswith("amr.max_grid_size_x")][0].split("=")[1].split()
    assert int(mgs[-1]) == p.fire_box and len(mgs) == len(p.levels)
    off = _plan(phenomenon="downslope_windstorm")
    assert not off.fire and not any("erf.fire" in ln for ln in writers.inputs_deck(off, "t"))
    on = _plan(phenomenon="downslope_windstorm", fire=True)
    assert on.fire and on.fire_ratio == 10
    assert not _plan(fire=False).fire


# ---- terrain ------------------------------------------------------------------------------

def test_gaussian_hill_slope():
    # h = H exp(-r^2 / (2 s^2)) has its steepest slope H / s exp(-1/2) at r = s.
    H, s = 400.0, 1500.0
    lat0, lon0 = 40.0, -105.0
    cs = 0.0005
    n = 241
    mlat = 111195.0
    mlon = mlat * math.cos(math.radians(lat0))
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "hill.asc")
        with open(path, "w") as f:
            f.write("ncols %d\nnrows %d\nxllcenter %r\nyllcenter %r\ncellsize %r\nNODATA_value -9999\n"
                    % (n, n, lon0 - cs * (n // 2), lat0 - cs * (n // 2), cs))
            for jj in range(n - 1, -1, -1):
                dy = (jj - n // 2) * cs * mlat
                vals = [H * math.exp(-(((ii - n // 2) * cs * mlon) ** 2 + dy ** 2) / (2 * s * s))
                        for ii in range(n)]
                f.write(" ".join("%.4f" % v for v in vals) + "\n")
        dem = AsciiGridDEM(path)
        p = build(Settings(lat=lat0, lon=lon0, phenomenon="hill_flow", target_dx=100.0,
                           max_base_dx=900.0, outer=(20.0, 20.0), finest=(4.0, 4.0), fetch=0.0))
        lv = p.levels[-1]
        h, water, dxs = level_terrain(p.frame, lv, dem, lat0)
        assert dxs == lv.dx
        st = slope_stats(h, dxs, 0.3, 0.6, water)
    _close(st["max_slope"], H / s * math.exp(-0.5), 0.02 * H / s, "Gaussian hill max slope")
    _close(st["relief"] + st["h_min"], H, 0.01 * H, "hill top")


def test_slope_warnings():
    # a ramp of slope 0.8 everywhere: every level above the 0.6 limit is warned about
    def ramp(lat, lon):
        return 0.8 * (lat - 39.0) * 111195.0
    p = build(Settings(lat=40.0, lon=-105.0, phenomenon="hill_flow", target_dx=300.0,
                       max_base_dx=900.0, outer=(30.0, 30.0), finest=(5.0, 5.0), fetch=0.0))
    stats, warns = check(p, AnalyticDEM(ramp, res=5.0), 0.3, 0.6, 150.0)
    for st in stats:
        _close(st["max_slope"], 0.8, 0.01, "ramp slope at level %d" % st["lev"])
    assert sum("steepest slope is 0.80" in w and "may fail" in w for w in warns) == len(p.levels), warns
    stats, warns = check(p, AnalyticDEM(lambda la, lo: 0.0 * la + 100.0), 0.3, 0.6, 150.0)
    assert not any("steepest" in w for w in warns), warns


def test_despike_and_mismatch():
    import numpy as np
    h = np.full((9, 9), 100.0)
    h[4, 4] = 5000.0
    out, n = despike(h, 200.0)
    assert n == 1 and out[4, 4] == 100.0, (n, out[4, 4])
    out, n = despike(h + np.arange(9)[None, :] * 50.0, 200.0)   # a smooth ramp is kept
    assert n == 1
    # one 500 m cell in the relaxation zone (one row in, so the window stays inside), 3 x 3
    # window at dx 1 km: 500 - 500 / 9
    h0 = np.zeros((20, 20))
    h0[1, 10] = 500.0
    _close(boundary_mismatch(h0, 1000.0, 5), 500.0 - 500.0 / 9.0, 1e-9, "boundary mismatch")
    h0[:] = 0.0
    h0[10, 10] = 500.0                                         # inside, outside the zone
    _close(boundary_mismatch(h0, 1000.0, 5), 0.0, 1e-12, "mismatch away from the edge")


def _png(rows, nch, filters):
    """A test encoder: 8-bit PNG with one filter type per row."""
    width = len(rows[0]) // nch
    raw = bytearray()
    prev = bytearray(len(rows[0]))
    for cur, f in zip(rows, filters):
        out = bytearray()
        for c, x in enumerate(cur):
            a = cur[c - nch] if c >= nch else 0
            b = prev[c]
            cc = prev[c - nch] if c >= nch else 0
            if f == 0:
                pred = 0
            elif f == 1:
                pred = a
            elif f == 2:
                pred = b
            elif f == 3:
                pred = (a + b) >> 1
            else:
                pa, pb, pc = abs(b - cc), abs(a - cc), abs(a + b - 2 * cc)
                pred = a if (pa <= pb and pa <= pc) else (b if pb <= pc else cc)
            out.append((x - pred) & 255)
        raw += bytes([f]) + out
        prev = cur

    def chunk(t, d):
        return struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d) & 0xffffffff)
    color = {1: 0, 3: 2, 4: 6}[nch]
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, len(rows), 8, color, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(bytes(raw))) + chunk(b"IEND", b""))


def test_png_decoder():
    rnd = random.Random(7)
    for nch in (1, 3, 4):
        rows = [bytearray(rnd.randrange(256) for _ in range(11 * nch)) for _ in range(10)]
        filters = [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]
        img = decode_png(_png(rows, nch, filters))
        assert img.shape == (10, 11, nch), img.shape
        for r in range(10):
            assert bytes(img[r].ravel()) == bytes(rows[r]), "row %d (filter %d, %d channels)" % (
                r, filters[r], nch)
    # Paeth ties: left a = 14, up b = 2, up-left c = 10 give pb = pc = 4 < pa = 8, so the
    # predictor is b (PNG specification: b wins a tie with c)
    rows = [bytearray([10, 2]), bytearray([14, 77])]
    img = decode_png(_png(rows, 1, [0, 4]))
    assert bytes(img[1].ravel()) == bytes(rows[1]), "Paeth tie b = c: %s" % list(img[1].ravel())


# ---- the deck -----------------------------------------------------------------------------

def test_deck_passes_checker():
    """The written deck has no checker ERROR or WARNING apart from the met_em files it names."""
    p = _plan(start="2025-01-07_12")
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "inputs")
        with open(path, "w") as f:
            f.write("\n".join(writers.inputs_deck(p, "selftest")) + "\n")
        res = subprocess.run([sys.executable, os.path.join(TOOLS, "check_erf_inputs.py"), path],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    bad = [ln for ln in res.stdout.splitlines()
           if (ln.startswith("[ERROR]") and "nc_init_file" not in ln) or ln.startswith("[WARNING]")]
    assert not bad, "checker findings:\n" + res.stdout
    assert "[ERROR] erf.nc_init_file_3" in res.stdout, "the checker did not look at level 3's file"


def _deck_value(deck, key):
    line = [ln for ln in deck if ln.split("=")[0].strip() == key][0]
    return line.split("=")[1].split("#")[0].split()


def test_deck_cell_size_exact():
    # ERF asserts |prob_hi / n_cell - met_em DX| < 1e-3 m; a 603.2 m base once lost digits
    p = build(Settings(lat=40.0, lon=-100.0, phenomenon="wind_farm", target_dx=37.7, ratio=2))
    deck = writers.inputs_deck(p, "t")
    hi = [float(v) for v in _deck_value(deck, "geometry.prob_hi")]
    n = [int(v) for v in _deck_value(deck, "amr.n_cell")]
    nl = writers.namelist_wps(p)
    dx = float([ln for ln in nl if ln.strip().startswith("dx =")][0].split("=")[1].strip(" ,"))
    for d in range(2):
        _close(hi[d] / n[d], dx, 1e-6, "prob_hi / n_cell against the namelist dx")
    _close(dx, p.base_dx, 1e-6, "namelist dx")


def test_deck_terrain_smoothing():
    # ERF aborts with terrain_smoothing != 0 on a refined level of a Metgrid run
    for kw in ({}, dict(phenomenon="hurricane_landfall")):
        p = _plan(**kw)
        deck = writers.inputs_deck(p, "t")
        assert _deck_value(deck, "erf.terrain_smoothing") == ["0"], kw
        assert not any(ln.startswith("prob.") or "real_set_width" in ln for ln in deck)


def test_time_step_erf_limit():
    # dt <= 0.8 dx / (m (sqrt(2) umax + 2 c / 3)), at most 6 s per km, the largest NICE_DT
    for dx, umax, m in ((900.0, 40.0, 1.0), (1000.0, 70.0, 1.02), (300.0, 20.0, 1.0), (100.0, 40.0, 1.0)):
        lim = min(6e-3 * dx, ERF_FIXED_DT_CFL * dx / (m * (math.sqrt(2) * umax + 2 * SOUND_SPEED / 3)))
        dt = base_time_step(dx, umax, m)
        assert dt <= lim + 1e-12, (dx, umax, dt, lim)
        bigger = [d for d in NICE_DT if d > dt]
        assert not bigger or bigger[0] > lim, (dx, umax, dt, lim)
    assert base_time_step(900.0, 40.0) == 2.0     # 0.8 * 900 / (40 sqrt(2) + 2 * 370 / 3) = 2.37
    assert all(abs(900.0 / d - round(900.0 / d)) < 1e-9 for d in NICE_DT), "15-minute output"


def test_vertical_edge_cases():
    nz, q, dz_top = vertical_grid(40000.0, 0.5, 200.0, 1.005)       # more than 1024 levels
    assert nz % 8 == 0 and nz > 1024 and 1.0 < q <= 1.005
    depth = 0.5 * (q ** nz - 1) / (q - 1)
    _close(depth, 40000.0, 1e-6 * 40000.0, "depth of the long grid")
    nz, q, dz = vertical_grid(1000.0, 150.0, 300.0, 1.08)           # uniform
    assert nz == 8 and q == 1.0, (nz, q)
    _close(dz, 125.0, 1e-12, "uniform dz = top / nz")
    p = _plan(top=1.0, dz0=150.0)
    assert p.dz0 == 125.0 and "erf.initial_dz" not in " ".join(writers.inputs_deck(p, "t"))


def test_vertical_grid_erf_newton():
    # ERF's geometric_stretch: Newton from 1.03, at most 50 steps, |F| <= 1e-8 m.  From 1.03 a
    # grid with q near 1.16 (48 cells of 0.5 m to 4 km) takes more than 50 steps; the advisor
    # adds cells until it converges
    def erf_converges(dz0, top, nz):
        s = 1.03
        f = dz0 * (s ** nz - 1) / (s - 1) - top
        it = 0
        while abs(f) > 1e-8 and it < 50:
            s -= f / (dz0 * (nz * s ** (nz - 1) * (s - 1) - s ** nz + 1) / (s - 1) ** 2)
            s = max(s, 1 + 1e-8)
            f = dz0 * (s ** nz - 1) / (s - 1) - top
            it += 1
        return abs(f) <= 1e-8
    assert not erf_converges(0.5, 4000.0, 48)
    nz, q, _ = vertical_grid(4000.0, 0.5, 1000.0, 1.2)
    assert nz % 8 == 0 and nz > 48 and erf_converges(0.5, 4000.0, nz), (nz, q)
    for top, dz0, ms in ((8000.0, 0.5, 1.3), (8000.0, 1.0, 1.2), (20000.0, 2.0, 1.3)):
        nz, q, _ = vertical_grid(top, dz0, 2000.0, ms)
        assert erf_converges(dz0, top, nz), (top, dz0, ms, nz)
    # uniform grid: top / nz must also meet the top-cell limit
    nz, q, dz = vertical_grid(1000.0, 150.0, 100.0, 1.08)
    assert q == 1.0 and dz <= 100.0 and nz == 16, (nz, q, dz)


def test_fire_ratio_two():
    # the fire level may grow past its parent before the parents are padded
    for kw in (dict(ratio=2, wind_from=0.0, lat=38.5, lon=-97.5),
               dict(ratio=2, fire_dx=3.0, target_dx=10.0, lat=39.0, lon=-105.0)):
        p = _plan(**kw)
        fl = p.levels[-1]
        assert fl.nx % p.fire_box == 0 and fl.ny % p.fire_box == 0
        assert fl.margin >= p.settings.nest_margin, (kw, fl.margin)


def test_cell_slope_one_cell():
    import numpy as np
    # a 100 m step between two cells 10 m apart: the grid sees 10, central differences 5
    h = np.array([[0.0, 0.0, 100.0, 100.0]] * 3)
    s = cell_slope(h, 10.0)
    _close(float(s.max()), 10.0, 1e-12, "one-cell slope of a step")
    _close(float(s[1, 0]), 0.0, 1e-12, "flat cell")
    _close(float(cell_slope(h.T, 10.0).max()), 10.0, 1e-12, "one-cell slope of a step in y")
    # a plane of slope 0.3 in x and 0.4 in y: 0.5 everywhere
    y, x = np.mgrid[0:5, 0:6] * 20.0
    _close(float(cell_slope(0.3 * x + 0.4 * y, 20.0).min()), 0.5, 1e-12, "plane slope")


def test_slope_factor_warning():
    # a 500 m band of slope 0.2 rising 100 m north of the point, on hill_flow (dz0 3 m): over the
    # top of the band the first cell is squeezed to dz0 (top - 100) / top, and
    # alpha = 0.2 dx / that > 1 on the 33 m level
    def ramp(lat, lon):
        import numpy as np
        return 0.2 * np.clip((lat - 40.0) * 111195.0, 0.0, 500.0)
    p = build(Settings(lat=40.0, lon=-105.0, phenomenon="hill_flow", outer=(30.0, 30.0),
                       finest=(2.0, 2.0), fetch=0.0))
    stats, warns = check(p, AnalyticDEM(ramp, res=5.0), 0.3, 0.6, 150.0)
    fl = p.levels[-1]
    _close(stats[-1]["h_max"], 100.0, 2.0, "highest ground")
    dz1 = p.dz0 * (p.top - stats[-1]["h_max"]) / p.top
    _close(stats[-1]["alpha"], 0.2 * fl.dx / dz1, 0.02, "slope factor")
    assert stats[-1]["alpha"] > 0.2 * fl.dx / p.dz0 * (1 + 1e-6), "first cell not squeezed"
    assert any(w.startswith("Level %d: the terrain slope factor" % fl.lev) for w in warns), warns
    assert not any("leaves" in w for w in warns), warns
    # a 5 km plateau under a top of a few km leaves no room for the boundary layer and damping
    p = build(Settings(lat=40.0, lon=-105.0, phenomenon="hill_flow", outer=(30.0, 30.0),
                       finest=(2.0, 2.0), fetch=0.0, top=6.0))
    stats, warns = check(p, AnalyticDEM(lambda lat, lon: 5000.0 + 0.0 * lat, res=5.0), 0.3, 0.6, 150.0)
    assert any("leaves" in w and "raise --top" in w for w in warns), warns


def test_damaged_tiles():
    good = _png([bytearray([128, 0, 0] * 4)] * 4, 3, [0] * 4)
    for bad in (good[:45], good[:33] + b"\x00" * 40):   # cut inside IDAT; no IDAT
        try:
            decode_png(bad)
        except TerrainError:
            pass
        else:
            raise AssertionError("damaged PNG accepted")
    with tempfile.TemporaryDirectory() as d:
        dem = TerrariumDEM(cache_dir=d, fetch=lambda url: good[:45])
        try:
            dem._tile(3, 1, 2)
        except TerrainError:
            pass
        else:
            raise AssertionError("damaged download accepted")
        assert not os.path.exists(os.path.join(d, "3", "1", "2.png")), "damaged tile was cached"
        dem = TerrariumDEM(cache_dir=d, fetch=lambda url: good)
        _close(float(dem._tile(3, 1, 2)[0, 0]), 0.0, 1e-9, "terrarium 128,0,0 is 0 m")
        assert os.path.exists(os.path.join(d, "3", "1", "2.png"))


def test_outside_hrrr_refused():
    # inside HRRR, but the 150 km outer domain reaches past its south-west corner
    _raises(lambda: _plan(lat=21.3, lon=-122.5), "outside the HRRR grid")


def test_terrain_tiles_antimeridian():
    import numpy as np
    good = _png([bytearray([128, 0, 0] * 4)] * 4, 3, [0] * 4)
    urls = []
    with tempfile.TemporaryDirectory() as d:
        dem = TerrariumDEM(cache_dir=d, fetch=lambda url: urls.append(url) or good)
        dem.prepare(np.array([52.0, 52.1]), np.array([179.9, -179.9]), 1000.0)
    assert dem.zoom >= 8 and len(urls) <= 4, "zoom %d, %d tiles fetched" % (dem.zoom, len(urls))


def test_boundary_mismatch_sampled():
    import numpy as np
    from . import terrain as T
    p = _plan()
    dem = AnalyticDEM(lambda la, lo: 500.0 * (np.asarray(lo) > -118.6) * (np.asarray(la) < 34.3))
    keep = {}
    old = T.MAX_SAMPLES
    T.MAX_SAMPLES = 2000            # level 0 (168 x 168) is sampled every 4 cells
    try:
        st, _ = T.check(p, dem, keep=keep)
    finally:
        T.MAX_SAMPLES = old
    h, _, dxs = keep[0]
    assert dxs > p.base_dx
    rim = max(1, int(math.ceil(p.settings.real_width * p.base_dx / dxs)))
    _close(st[0]["boundary_mismatch"], boundary_mismatch(h, dxs, rim), 1e-9, "mismatch on sampled spacing")


def test_even_box_multiple_of_ratio():
    for n, r in ((123, 3), (86, 2), (720, 3), (97, 3)):
        b = writers._even_box(n, r)
        assert b % r == 0 and 24 <= b <= 128, (n, r, b)
    assert writers._even_box(168, 1) == 56


def test_outdir_is_a_file():
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "plain_file")
        open(f, "w").close()
        r = subprocess.run([sys.executable, _tool(), "--lat", "40", "--lon", "-105", "--phenomenon", "hill_flow",
                            "--dem", "none", "--outdir", f], capture_output=True, text=True)
    assert r.returncode == 2 and "Traceback" not in r.stderr and "--outdir" in r.stderr, r.stderr[-300:]


def _tool():
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "erfdomainadvisor.py")


def test_fire_flags_exclusive():
    r = subprocess.run([sys.executable, _tool(), "--lat", "40", "--lon", "-105", "--phenomenon",
                        "wildfire", "--dem", "none", "--fire", "--no-fire"], capture_output=True, text=True)
    assert r.returncode == 2 and "not allowed with" in r.stderr, r.stderr[-300:]


def test_start_whole_hour():
    assert writers.parse_start("2025-01-07_12").hour == 12
    assert writers.parse_start("2025-01-07_12:00:00").hour == 12
    for bad in ("2025-01-07", "2025-01-07_12:30:00", "2025-01-07 12:00:05", "07/01/2025"):
        try:
            writers.parse_start(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("--start %s accepted" % bad)


def test_shallow_top_warned():
    # convective_bl_les: a 1.5 km top under a deep boundary layer leaves no room for damping
    kw = dict(phenomenon="convective_bl_les", wind_from=None, lat=36.6, lon=-97.5)
    p = _plan(**kw)
    assert not any("model top" in w for w in p.warnings), p.warnings
    p = _plan(top=(p.zi + 100.0) / 1000.0, **kw)
    assert any("model top" in w and "Raise --top" in w for w in p.warnings), p.warnings


def test_gray_zone_text():
    # the warnings state the same bounds as regime(): dx = 0.2 zi is in the gray zone
    p = _plan(phenomenon="convective_bl_les", wind_from=None, lat=36.6, lon=-97.5)
    gray = [w for w in p.warnings if "gray zone" in w]
    assert gray and all("0.2 zi <= dx < zi" in w for w in gray), gray


TESTS = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
# tests that need numpy; the others run without it
NEEDS_NUMPY = {"test_gaussian_hill_slope", "test_slope_warnings", "test_despike_and_mismatch",
               "test_png_decoder", "test_cell_slope_one_cell", "test_slope_factor_warning",
               "test_damaged_tiles", "test_terrain_tiles_antimeridian", "test_boundary_mismatch_sampled"}



def main():
    try:
        import numpy  # noqa: F401
        tests = TESTS
    except ImportError:
        tests = [t for t in TESTS if t.__name__ not in NEEDS_NUMPY]
        print("self-test: numpy is not installed; the %d terrain tests are skipped"
              % (len(TESTS) - len(tests)))
    failed = 0
    for t in tests:
        try:
            t()
            print("PASS %s" % t.__name__)
        except Exception as exc:      # report every failure, then fail the run
            failed += 1
            print("FAIL %s: %s: %s" % (t.__name__, type(exc).__name__, exc))
    print("self-test: %s" % ("%d tests passed" % len(tests) if not failed
                             else "%d of %d tests FAILED" % (failed, len(tests))))
    return 1 if failed else 0
