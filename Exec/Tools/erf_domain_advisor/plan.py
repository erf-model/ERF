"""Domain, refinement levels, vertical grid, time step and cost for one point and phenomenon.

The rules, all best guesses (see README.md):

* Base cell size: the finest cell size times ratio**L, with L the largest number of levels that
  keeps the base cell at or below --max-base-dx (default 1000 m, a third of the HRRR spacing).
* Each refined box sits at least --nest-margin parent cells (default 12) inside its parent.
* The outer domain keeps the first refined box at least erf.real_width base cells plus
  --hrrr-buffer HRRR cells (default 5, 15 km) away from its edges, and is at least the
  phenomenon's outer size.
* The finest box gets the phenomenon's upwind fetch: all of it upwind with --wind-from, half of
  it on every side without.
* Boxes are snapped outward onto their parent's cells, so each nest is a WPS nest with
  e_we - 1 and e_sn - 1 divisible by the ratio; ERF reads its place from the met_em file.
* The vertical grid has a first cell of dz0 and nz cells; ERF solves the stretching ratio
  (erf.initial_dz without erf.grid_stretching_ratio).  nz is the smallest multiple of 8 that
  keeps the ratio at or below --max-stretch and the top cell at or below the phenomenon's dz_top;
  when nz cells of dz0 already reach the top, the grid is uniform.
* The base time step is the largest step that divides 15 minutes and stays at or below both
  6 s per km and the Courant number of 0.8 that ERF enforces on a fixed dt (sound speed over 3
  per direction plus the phenomenon's largest wind, times the largest map factor).
* With fire, erf.fire.grid_ratio C is the smallest integer that makes the fire cell (finest dx
  over C) at most --fire-dx (default 10 m).  ERF needs every x,y box length on the fire level
  to be a multiple of C, so the finest box is grown to a multiple of a box length B (a multiple
  of C and of the ratio, near 64 cells) and the deck sets amr.max_grid_size_x/y = B there.
"""

import math
from dataclasses import dataclass, field
from typing import List, Optional

from .phenomena import PHENOMENA
from .projection import (HRRR_DX, DomainFrame, Lambert, WPS_RADIUS, choose_projection,
                         hrrr_margin_cells)

# Time steps [s] that divide 15 minutes, so hourly and 15-minute output fall on a step.
NICE_DT = [0.25, 0.5, 1.0, 1.5, 2.0, 2.4, 2.5, 3.0, 3.6, 4.0, 4.5, 5.0, 6.0, 7.2, 7.5, 9.0,
           10.0, 12.0, 15.0, 18.0, 20.0]
SOUND_SPEED = 370.0       # [m/s], sqrt(1.4 R T) at about 340 K: margin over near-surface air; hot fire cells are faster
ERF_FIXED_DT_CFL = 0.8    # ERF lowers a fixed dt above this Courant number (default erf.cfl)


class PlanError(Exception):
    """An input outside its range, or a plan that cannot be built."""


@dataclass
class Settings:
    lat: float
    lon: float
    phenomenon: str
    target_dx: Optional[float] = None
    ratio: int = 3
    max_base_dx: float = 1000.0
    outer: Optional[tuple] = None
    finest: Optional[tuple] = None
    fetch: Optional[float] = None
    top: Optional[float] = None
    dz0: Optional[float] = None
    hours: Optional[float] = None
    wind_from: Optional[float] = None
    real_width: int = 5
    nest_margin: int = 12
    hrrr_buffer: int = 5
    max_stretch: float = 1.08
    cost_per_cell_step: float = 2.0e-6
    bytes_per_cell: float = 2000.0
    start: Optional[str] = None
    fire: Optional[bool] = None      # None: the phenomenon's default
    fire_dx: float = 10.0


@dataclass
class Level:
    lev: int
    dx: float
    x0: float          # box in ERF level-0 coordinates [m]
    y0: float
    x1: float
    y1: float
    nx: int            # cells at this level
    ny: int
    i_start: int = 0   # 0-based start in the parent's index space (levels > 0)
    j_start: int = 0
    margin: int = 0    # smallest distance to the parent's edge, parent cells (levels > 0)
    regime: str = ""   # "PBL scheme", "gray zone" or "LES"
    dt: float = 0.0


@dataclass
class Plan:
    settings: Settings
    title: str
    notes: str
    moisture: str
    zi: float
    umax: float
    target_dx: float
    base_dx: float
    ratio: int
    levels: List[Level]
    frame: DomainFrame
    hrrr_projection: bool
    top: float
    dz0: float
    nz: int
    stretch: float
    dz_top: float
    hours: float
    fetch: float
    point_xy: tuple
    hrrr_margin: float
    warnings: List[str] = field(default_factory=list)
    cells: List[int] = field(default_factory=list)
    core_hours: float = 0.0
    memory_gb: float = 0.0
    fire: bool = False
    fire_ratio: int = 0    # erf.fire.grid_ratio
    fire_box: int = 0      # x,y box length on the fire level, a multiple of fire_ratio


def _check(cond, msg):
    if not cond:
        raise PlanError(msg)


def validate(s: Settings):
    """Range checks; a value out of range stops with a message naming the option (no clamping)."""
    _check(-80.0 <= s.lat <= 80.0, "--lat must lie in [-80, 80] (a Lambert conformal grid), got %g" % s.lat)
    _check(-180.0 <= s.lon <= 360.0, "--lon must lie in [-180, 360], got %g" % s.lon)
    _check(s.phenomenon in PHENOMENA, "--phenomenon %r is not known; see --list" % s.phenomenon)
    if s.target_dx is not None:
        _check(1.0 <= s.target_dx <= 3000.0, "--target-dx must lie in [1, 3000] m, got %g" % s.target_dx)
    _check(s.ratio in (2, 3), "--ratio must be 2 or 3, got %s" % s.ratio)
    _check(300.0 <= s.max_base_dx <= 3000.0,
           "--max-base-dx must lie in [300, 3000] m, got %g" % s.max_base_dx)
    for name, v in (("--outer", s.outer), ("--finest", s.finest)):
        if v is not None:
            _check(len(v) == 2 and all(1.0 <= a <= 5000.0 for a in v),
                   "%s needs two lengths in [1, 5000] km, got %s" % (name, v))
    if s.fetch is not None:
        _check(0.0 <= s.fetch <= 200.0, "--fetch must lie in [0, 200] km, got %g" % s.fetch)
    if s.top is not None:
        _check(1.0 <= s.top <= 40.0, "--top must lie in [1, 40] km, got %g" % s.top)
    if s.dz0 is not None:
        _check(0.5 <= s.dz0 <= 200.0, "--dz0 must lie in [0.5, 200] m, got %g" % s.dz0)
    if s.hours is not None:
        _check(0.25 <= s.hours <= 240.0, "--hours must lie in [0.25, 240], got %g" % s.hours)
    if s.wind_from is not None:
        _check(0.0 <= s.wind_from < 360.0, "--wind-from must lie in [0, 360) degrees, got %g" % s.wind_from)
    _check(1 <= s.real_width <= 50, "--real-width must lie in [1, 50] cells, got %d" % s.real_width)
    _check(4 <= s.nest_margin <= 100, "--nest-margin must lie in [4, 100] cells, got %d" % s.nest_margin)
    _check(0 <= s.hrrr_buffer <= 100, "--hrrr-buffer must lie in [0, 100] cells, got %d" % s.hrrr_buffer)
    _check(1.0 < s.max_stretch <= 1.3, "--max-stretch must lie in (1, 1.3], got %g" % s.max_stretch)
    _check(0.5 <= s.fire_dx <= 100.0, "--fire-dx must lie in [0.5, 100] m, got %g" % s.fire_dx)
    _check(1e-9 <= s.cost_per_cell_step <= 1e-3,
           "--cost-per-cell-step must lie in [1e-9, 1e-3] core-s, got %g" % s.cost_per_cell_step)


def damping_depth(top):
    """Depth [m] of the Rayleigh damping layer below the model top that the deck writes."""
    return min(0.3 * top, 8000.0)


def level_chain(target_dx, ratio, max_base_dx):
    """(number of refined levels, base cell size)."""
    if target_dx >= max_base_dx:
        return 0, target_dx
    nlev = int(math.floor(math.log(max_base_dx / target_dx) / math.log(ratio) + 1e-9))
    if nlev > 6:
        raise PlanError("--target-dx %g m needs %d refined levels from a %g m base; at most 6 are "
                        "allowed" % (target_dx, nlev, max_base_dx))
    return nlev, target_dx * ratio ** nlev


def vertical_grid(top, dz0, dz_top, max_stretch):
    """(nz, ratio, top cell height) of a geometric grid of nz cells that fills top from dz0.
    nz is a multiple of 8 (it splits evenly into GPU-friendly boxes); when nz uniform cells of
    dz0 already reach the top, the grid is uniform with dz = top / nz, as ERF builds it."""
    if dz0 * 4 > top:
        raise PlanError("--dz0 %g m is too thick for a %g m deep domain" % (dz0, top))
    nz = 4
    while nz < 2000:
        if nz * dz0 >= top and top / nz <= dz_top:
            break
        q = solve_stretch(nz, dz0, top)
        if q <= max_stretch and dz0 * q ** (nz - 1) <= dz_top:
            break
        nz += 1
    else:
        raise PlanError("no vertical grid of fewer than 2000 cells meets dz0 = %g m and a top cell "
                        "of at most %g m; raise --dz0 or --max-stretch, or lower --top"
                        % (dz0, dz_top))
    # more cells only lower q and the top cell
    nz = 8 * int(math.ceil(nz / 8.0))
    if nz * dz0 >= top:
        return nz, 1.0, top / nz
    # ERF finds q by Newton's method from 1.03 in at most 50 steps and stops if it has not
    # converged; from 1.03 a large q needs many steps, and more cells lower q
    while not erf_newton_converges(dz0, top, nz):
        nz += 8
        if nz >= 2000:
            raise PlanError("ERF's stretching solve does not converge for dz0 = %g m and a %g m "
                            "deep domain with fewer than 2000 cells; raise --dz0" % (dz0, top))
    q = solve_stretch(nz, dz0, top)
    return nz, q, dz0 * q ** (nz - 1)


ERF_NEWTON_TOL = 1.0e-8     # erf_init_zlevels passes 1e-8 m in double precision
ERF_NEWTON_STEPS = 45       # ERF allows 50; 5 kept in hand for round-off


def erf_newton_converges(dz0, top, nz, tol=ERF_NEWTON_TOL, steps=ERF_NEWTON_STEPS):
    """True when erf_grid_utils::geometric_stretch (Source/Utils/ERF_GridUtils.H) converges for
    this grid within steps iterations: the same Newton iteration from SFact = 1.03 on
    F = dz0 (S**nz - 1) / (S - 1) - top, with S kept at or above 1 + tol."""
    if nz * dz0 >= top:
        return True
    s = 1.03
    try:
        f = dz0 * (s ** nz - 1.0) / (s - 1.0) - top
        for _ in range(steps):
            if abs(f) <= tol:
                return True
            df = dz0 * (nz * s ** (nz - 1) * (s - 1.0) - s ** nz + 1.0) / (s - 1.0) ** 2
            s = max(s - f / df, 1.0 + tol)
            f = dz0 * (s ** nz - 1.0) / (s - 1.0) - top
    except OverflowError:
        return False
    return abs(f) <= tol


def solve_stretch(nz, dz0, top):
    """Ratio q > 1 with dz0 (q**nz - 1) / (q - 1) = top, for nz dz0 < top.  Bisection on ln q,
    with the sum written as expm1(nz ln q) / expm1(ln q) so that large nz cannot overflow."""
    def depth(lq):
        return dz0 * math.expm1(nz * lq) / math.expm1(lq)
    lo, hi = 1e-12, max(1e-6, math.log(top / dz0) / max(nz - 1, 1))
    while depth(hi) < top:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if depth(mid) < top:
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (lo + hi))


def fire_grid(finest_dx, fire_dx, ratio):
    """(erf.fire.grid_ratio C, box length B) for the fire level: C is the smallest integer with
    finest_dx / C <= fire_dx; B is a multiple of C and of ratio (the level's x,y lengths are
    multiples of ratio) near 64 cells."""
    c = max(1, int(math.ceil(finest_dx / fire_dx - 1e-9)))
    if c > 16:
        raise PlanError("fire cells of at most --fire-dx %g m need erf.fire.grid_ratio %d on a %g m "
                        "finest level; at most 16 is allowed. Lower --target-dx or raise --fire-dx."
                        % (fire_dx, c, finest_dx))
    lcm = c * ratio // math.gcd(c, ratio)
    return c, lcm * max(1, int(round(64.0 / lcm)))


def _grow_to_multiple(i0, i1, m):
    """Grow [i0, i1) by as few cells as possible, about equally on both sides, so that its
    length is a multiple of m."""
    extra = (-(i1 - i0)) % m
    return i0 - extra // 2, i1 + extra - extra // 2


def regime(dx, zi):
    """Turbulence regime of a cell size against the boundary-layer depth zi (best guess)."""
    if dx >= zi:
        return "PBL scheme"
    if dx < 0.2 * zi:
        return "LES"
    return "gray zone"


def base_time_step(base_dx, umax, map_factor=1.0):
    """Largest NICE_DT at or below min(6 s per km, ERF's fixed-dt limit).

    ERF lowers a fixed dt whose Courant number exceeds min(erf.cfl, 0.8) (fixed_dt_cfl_limit,
    on by default).  With acoustic substeps its rate is (|u| + c/3)/dx + (|v| + c/3)/dy, map
    factor m included; the worst case is u = v = umax/sqrt(2):
        dt <= 0.8 dx / (m (sqrt(2) umax + 2 c / 3)),   c = SOUND_SPEED.
    The 6 s per km cap is the usual WRF practice."""
    rate = map_factor * (math.sqrt(2.0) * umax + 2.0 * SOUND_SPEED / 3.0) / base_dx
    limit = min(6.0 * base_dx / 1000.0, ERF_FIXED_DT_CFL / rate)
    good = [d for d in NICE_DT if d <= limit + 1e-12]
    if not good:
        raise PlanError("the time step for dx = %g m would be below %g s" % (base_dx, NICE_DT[0]))
    return good[-1]


def _downstream_unit(proj, lat, lon, wind_from):
    """Unit vector, in grid x/y, of the direction the wind blows toward."""
    th = math.radians(wind_from)
    east, north = -math.sin(th), -math.cos(th)
    d = 1e-3
    x0, y0 = proj.forward(lat, lon)
    x1, y1 = proj.forward(lat + d * north, lon + d * east / max(math.cos(math.radians(lat)), 1e-6))
    vx, vy = x1 - x0, y1 - y0
    n = math.hypot(vx, vy)
    return vx / n, vy / n


def build(s: Settings) -> Plan:
    validate(s)
    p = PHENOMENA[s.phenomenon]
    lon = (s.lon + 180.0) % 360.0 - 180.0
    target = s.target_dx if s.target_dx is not None else p.target_dx
    outer = tuple(1000.0 * a for a in (s.outer or p.outer))
    finest = tuple(1000.0 * a for a in (s.finest or p.finest))
    fetch = 1000.0 * (s.fetch if s.fetch is not None else p.fetch)
    top = 1000.0 * (s.top if s.top is not None else p.top)
    dz0 = s.dz0 if s.dz0 is not None else p.dz0
    hours = s.hours if s.hours is not None else p.hours
    warnings = []

    nlev, base = level_chain(target, s.ratio, s.max_base_dx)
    if abs(base - round(base, -1)) <= 0.005 * base:
        # 33.3 m typed for 900/27 m: keep the base a round number of metres
        base = round(base, -1)
        target = base / s.ratio ** nlev
    if HRRR_DX / base > 5.0 + 1e-9:
        warnings.append("The base cell (%g m) is more than 5 times finer than HRRR (3 km); raise "
                        "--max-base-dx or choose a target that divides 900 m by powers of 3."
                        % base)
    dxs = [base / s.ratio ** l for l in range(nlev + 1)]
    fire = p.fire if s.fire is None else s.fire
    fire_ratio = fire_box = 0
    if fire:
        fire_ratio, fire_box = fire_grid(dxs[-1], s.fire_dx, s.ratio if nlev else 2)
        if dxs[-1] > 100.0:
            warnings.append("Fire on a %g m finest level: the plume and the fire-induced winds are "
                            "not resolved above about 100 m (a best guess); lower --target-dx."
                            % dxs[-1])

    truelat, stand_lon, in_hrrr = choose_projection(s.lat, lon)
    proj = Lambert(truelat, stand_lon, radius=WPS_RADIUS)
    if not in_hrrr:
        warnings.append("The point is outside the HRRR CONUS grid: HRRR cannot drive this case "
                        "(Alaska HRRR or GFS would be needed). The WPS domain uses the point's "
                        "own latitude and longitude for its Lambert projection.")
    if abs(s.lat) < 20.0:
        warnings.append("Near the equator a Mercator projection (map_proj = 'mercator') suits WPS "
                        "better than Lambert; change namelist.wps by hand.")

    # Region of interest (the finest box) around the point, in grid metres from the point.
    hx, hy = finest[0] / 2, finest[1] / 2
    box = [-hx, -hy, hx, hy]
    if s.wind_from is not None:
        ux, uy = _downstream_unit(proj, s.lat, lon, s.wind_from)
        if ux > 0:
            box[0] -= fetch * ux
        else:
            box[2] -= fetch * ux
        if uy > 0:
            box[1] -= fetch * uy
        else:
            box[3] -= fetch * uy
    else:
        ux = uy = 0.0
        box = [box[0] - fetch / 2, box[1] - fetch / 2, box[2] + fetch / 2, box[3] + fetch / 2]
        if fetch > 0 and nlev > 0:
            warnings.append("No --wind-from: the %g km fetch is split over all four sides of the "
                            "finest box. Give the wind direction to put all of it upwind." % (fetch / 1000))

    pad = 0.0
    for _attempt in range(8):
        # Boxes from the finest outward.
        # pad: extra room around the finest box for its parents, so that growing the fire level to
        # whole fire boxes leaves the nest margin intact (raised and rebuilt below when needed).
        seed = [box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad]
        boxes = {nlev: box}
        for l in range(nlev - 1, 0, -1):
            m = s.nest_margin * dxs[l]
            b = boxes[l + 1] if l < nlev - 1 else seed
            boxes[l] = [b[0] - m, b[1] - m, b[2] + m, b[3] + m]
        inner = boxes[1] if nlev >= 2 else seed
        buf = s.real_width * base + s.hrrr_buffer * HRRR_DX
        if nlev >= 1:
            buf += s.nest_margin * base
        need = [inner[0] - buf, inner[1] - buf, inner[2] + buf, inner[3] + buf]
        cx, cy = 0.5 * (box[0] + box[2]), 0.5 * (box[1] + box[3])
        shift = 0.15 * min(outer)
        cx, cy = cx - shift * ux, cy - shift * uy
        want = [cx - outer[0] / 2, cy - outer[1] / 2, cx + outer[0] / 2, cy + outer[1] / 2]
        ob = [min(need[0], want[0]), min(need[1], want[1]), max(need[2], want[2]), max(need[3], want[3])]
        m0 = fire_box if (fire and nlev == 0) else 2
        nx0 = m0 * int(math.ceil((ob[2] - ob[0]) / base / m0 - 1e-9))
        ny0 = m0 * int(math.ceil((ob[3] - ob[1]) / base / m0 - 1e-9))
        lx, ly = nx0 * base, ny0 * base
        ocx, ocy = 0.5 * (ob[0] + ob[2]), 0.5 * (ob[1] + ob[3])

        xp, yp = proj.forward(s.lat, lon)
        ref_lat, ref_lon = proj.inverse(xp + ocx, yp + ocy)
        frame = DomainFrame(ref_lat, ref_lon, lx, ly, truelat, stand_lon)
        px, py = frame.to_erf(s.lat, lon)
        # Grid metres from the point -> ERF coordinates.  The frame is centred on the outer centre,
        # so this offset is exact up to round-off, which is removed (to 1 mm) so that a box edge
        # on a cell face is not snapped differently on different levels.
        offx, offy = round(px, 3), round(py, 3)

        levels = [Level(0, base, 0.0, 0.0, lx, ly, nx0, ny0)]
        for l in range(1, nlev + 1):
            par = levels[-1]
            b = boxes[l]
            dxp = par.dx
            i0 = int(math.floor((b[0] + offx - par.x0) / dxp + 1e-9))
            j0 = int(math.floor((b[1] + offy - par.y0) / dxp + 1e-9))
            i1 = int(math.ceil((b[2] + offx - par.x0) / dxp - 1e-9))
            j1 = int(math.ceil((b[3] + offy - par.y0) / dxp - 1e-9))
            if fire and l == nlev:
                # fire level: x,y lengths multiples of the fire box length (in parent cells)
                i0, i1 = _grow_to_multiple(i0, i1, fire_box // s.ratio)
                j0, j1 = _grow_to_multiple(j0, j1, fire_box // s.ratio)
            margin = min(i0, j0, par.nx - i1, par.ny - j1)
            if margin < 1 and not (fire and l == nlev):
                # the fire level may grow past its parent; the retry below pads the parents
                raise PlanError("level %d does not fit inside level %d; enlarge --outer" % (l, l - 1))
            levels.append(Level(l, dxs[l], par.x0 + i0 * dxp, par.y0 + j0 * dxp,
                                par.x0 + i1 * dxp, par.y0 + j1 * dxp,
                                (i1 - i0) * s.ratio, (j1 - j0) * s.ratio, i0, j0, margin))

        need_margin = s.nest_margin if nlev >= 2 else 1
        if not fire or nlev == 0 or levels[nlev].margin >= need_margin:
            break
        pad += (need_margin - levels[nlev].margin) * dxs[nlev - 1]
    else:
        raise PlanError("could not fit the fire level's whole boxes inside level %d" % (nlev - 1))

    if fire:
        fl = levels[-1]
        ranks = (fl.nx // fire_box) * (fl.ny // fire_box)
        if ranks < 8:
            warnings.append("The fire level splits into %d boxes of %d x %d cells: the run can use "
                            "at most %d ranks before AMReX cuts the boxes and breaks the fire grid."
                            % (ranks, fire_box, fire_box, ranks))
    if nlev > 0 and levels[-1].dx < 900.0:
        warnings.append("WPS interpolates its default topography (GMTED2010 30 arc-seconds, about "
                        "900 m) onto the nests: add a finer data set (for example 3DEP or SRTM 1 "
                        "arc-second) to geog_data_path and GEOGRID.TBL and set geog_data_res, or "
                        "the nests get smoother terrain than the terrain check uses.")
    for lv in levels:
        lv.regime = regime(lv.dx, p.zi)
    if levels[0].regime == "gray zone":
        warnings.append("Level 0 (%g m) lies in the gray zone for a boundary layer %g m deep "
                        "(0.2 zi <= dx < zi); it uses the PBL scheme, as is usual near 1 km."
                        % (base, p.zi))
    gray = [lv.lev for lv in levels[1:] if lv.regime == "gray zone"]
    if gray:
        warnings.append("Level(s) %s lie in the gray zone (0.2 zi <= dx < zi, zi = %g m): the "
                        "boundary layer is neither parameterised well nor resolved. Keep those "
                        "boxes small; they only carry the flow to the next level."
                        % (", ".join(map(str, gray)), p.zi))
    if levels[-1].regime != "LES" and target < p.zi:
        warnings.append("The finest level (%g m) does not resolve the boundary layer (zi = %g m)."
                        % (target, p.zi))

    nz, stretch, dz_top = vertical_grid(top, dz0, p.dz_top, s.max_stretch)
    if top < p.zi + damping_depth(top):
        warnings.append("The model top (%g m) is below the boundary layer (%g m) plus the damping "
                        "layer the deck writes (%g m): the damping reaches into the boundary "
                        "layer. Raise --top." % (top, p.zi, damping_depth(top)))
    if stretch == 1.0 and abs(dz_top - dz0) > 1e-9:
        warnings.append("%d cells of dz0 = %g m overfill the %g m top: the grid is uniform with "
                        "dz = %g m." % (nz, dz0, top, dz_top))
        dz0 = dz_top
    if levels[-1].dx < dz0:
        warnings.append("The finest dx (%g m) is below dz0 (%g m); consider a smaller --dz0."
                        % (levels[-1].dx, dz0))

    map_max = max(frame.proj.map_factor(frame.to_latlon(x, y)[0]) for x, y in _edge_points(lx, ly, 8))
    dt0 = base_time_step(base, p.umax, map_max)
    for lv in levels:
        lv.dt = dt0 / s.ratio ** lv.lev

    margin = min(hrrr_margin_cells(*frame.to_latlon(x, y))
                 for x, y in _edge_points(lx, ly, 40))
    if in_hrrr and margin < 0:
        raise PlanError("part of the outer domain lies outside the HRRR grid (%.1f HRRR cells "
                        "beyond its edge), where metgrid has no data: move the point inland of "
                        "the HRRR edge or reduce --outer" % -margin)
    if in_hrrr and margin < s.hrrr_buffer:
        warnings.append("The outer domain comes within %.1f HRRR cells of the HRRR grid edge "
                        "(%s wanted)%s." % (margin, s.hrrr_buffer,
                                            "; part of it is outside HRRR" if margin < 0 else ""))

    plan = Plan(s, p.title, p.notes, p.moisture, p.zi, p.umax, target, base, s.ratio, levels, frame,
                in_hrrr, top, dz0, nz, stretch, dz_top, hours, fetch, (px, py), margin, warnings,
                fire=fire, fire_ratio=fire_ratio, fire_box=fire_box)
    seconds = hours * 3600.0
    total = 0.0
    for lv in levels:
        cells = lv.nx * lv.ny * nz
        plan.cells.append(cells)
        total += cells * seconds / lv.dt
    plan.core_hours = total * s.cost_per_cell_step / 3600.0
    plan.memory_gb = sum(plan.cells) * s.bytes_per_cell / 1e9
    return plan


def _edge_points(lx, ly, n):
    for k in range(n + 1):
        t = k / n
        yield t * lx, 0.0
        yield t * lx, ly
        yield 0.0, t * ly
        yield lx, t * ly


def erf_nest_box(levels, ratio, lev):
    """The fine box ERF builds for level lev from the met_em attributes, in level-lev index
    space: ERF's Tagging reads (i_parent_start - 1, nx) from the file, forms
    Box(ratio*is, ratio*is + nx - 1), coarsens it by the ratio and shifts it by the parent's
    lower corner.  Returned as (ilo, jlo, ihi, jhi) at level lev."""
    lo_i = lo_j = 0
    for l in range(1, lev + 1):
        lv = levels[l]
        ci0 = (ratio * lv.i_start) // ratio + lo_i
        cj0 = (ratio * lv.j_start) // ratio + lo_j
        lo_i, lo_j = ci0 * ratio, cj0 * ratio
    lv = levels[lev]
    return lo_i, lo_j, lo_i + lv.nx - 1, lo_j + lv.ny - 1
