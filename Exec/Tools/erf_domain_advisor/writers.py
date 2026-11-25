"""Report, namelist.wps, ERF inputs file and plan.json for a plan."""

import datetime as _dt
import json
import math

from .plan import damping_depth


def _g(v):
    return "%g" % float("%.6g" % v)


def _exact(v):
    """Enough digits that ERF recomputes the same cell size (ERF checks the met_em DX to 1e-3 m)."""
    return "%.10g" % v


def times(plan):
    """met_em times (hourly), or None without --start."""
    s = plan.settings.start
    if not s:
        return None
    t0 = parse_start(s)
    n = int(math.ceil(plan.hours))
    return [t0 + _dt.timedelta(hours=h) for h in range(n + 1)]


def parse_start(s):
    """First time of the run; HRRR analyses are hourly, so the hour is needed and minutes and
    seconds must be 0."""
    for fmt in ("%Y-%m-%d_%H", "%Y-%m-%d_%H:%M:%S", "%Y-%m-%dT%H", "%Y-%m-%d %H:%M:%S"):
        try:
            t = _dt.datetime.strptime(s, fmt)
        except ValueError:
            continue
        if t.minute or t.second:
            raise ValueError("--start %r: HRRR files are hourly; use a whole hour" % s)
        return t
    raise ValueError("--start %r: use YYYY-MM-DD_HH (or YYYY-MM-DD_HH:00:00)" % s)


def wps_name(t):
    return t.strftime("%Y-%m-%d_%H:%M:%S")


# --------------------------------------------------------------------------------------------

def report(plan, tstats=None, twarn=(), dem_name=None):
    s, f = plan.settings, plan.frame
    out = []
    out.append("ERF domain advisor (best guesses; review before use)")
    out.append("=" * 72)
    out.append("Point            : %.4f, %.4f" % (s.lat, (s.lon + 180) % 360 - 180))
    out.append("Phenomenon       : %s (%s)" % (s.phenomenon, plan.title))
    out.append("Notes            : %s" % plan.notes)
    out.append("")
    out.append("Resolution       : finest %s m, base %s m, %d refined level(s), ratio %d"
               % (_g(plan.target_dx), _g(plan.base_dx), len(plan.levels) - 1, plan.ratio))
    out.append("Outer domain     : %.1f x %.1f km, centre %.4f, %.4f"
               % (f.lx / 1000, f.ly / 1000, f.ref_lat, f.ref_lon))
    out.append("Point in ERF     : x = %.0f m, y = %.0f m" % plan.point_xy)
    out.append("Projection       : Lambert, truelat %g, stand_lon %g%s"
               % (f.truelat, f.stand_lon, " (same as HRRR)" if plan.hrrr_projection else ""))
    out.append("HRRR edge margin : %.1f HRRR cells from the outer domain" % plan.hrrr_margin)
    out.append("Vertical         : top %.1f km, nz %d, dz0 %s m, ratio %.4f, top cell %.0f m"
               % (plan.top / 1000, plan.nz, _g(plan.dz0), plan.stretch, plan.dz_top))
    out.append("Run              : %s h, base time step %s s (subcycled by %d per level)"
               % (_g(plan.hours), _g(plan.levels[0].dt), plan.ratio))
    if plan.fire:
        out.append("Fire             : on level %d, erf.fire.grid_ratio %d (fire cells %s m), "
                   "x,y boxes of %d cells" % (len(plan.levels) - 1, plan.fire_ratio,
                                              _g(plan.levels[-1].dx / plan.fire_ratio), plan.fire_box))
    else:
        out.append("Fire             : off (--fire adds an erf.fire.* block)")
    out.append("")
    out.append("%-4s %-8s %-23s %-14s %-9s %-11s %-10s" % ("lev", "dx [m]", "box x, y [km]", "cells",
                                                          "margin", "regime", "dt [s]"))
    for lv, cells in zip(plan.levels, plan.cells):
        box = "%.1f-%.1f, %.1f-%.1f" % (lv.x0 / 1e3, lv.x1 / 1e3, lv.y0 / 1e3, lv.y1 / 1e3)
        out.append("%-4d %-8s %-23s %-14s %-9s %-11s %-10s" % (
            lv.lev, _g(lv.dx), box, "%dx%dx%d" % (lv.nx, lv.ny, plan.nz),
            "-" if lv.lev == 0 else "%d" % lv.margin, lv.regime, _g(lv.dt)))
    out.append("")
    out.append("Cost (rough)     : %.3g cells, %.3g core-hours, %.3g GB"
               % (sum(plan.cells), plan.core_hours, plan.memory_gb))
    out.append("                   (%g core-s per cell-step and %g bytes per cell; both guesses)"
               % (s.cost_per_cell_step, s.bytes_per_cell))
    if tstats:
        out.append("")
        out.append("Terrain (%s, about %.0f m data)" % (dem_name, tstats[-1]["dem_res"]))
        out.append("%-4s %-8s %-10s %-10s %-10s %-10s %-8s" % ("lev", "dx [m]", "max slope",
                                                              "p99 slope", "relief [m]", "water", ""))
        for st in tstats:
            extra = ""
            if "boundary_mismatch" in st:
                extra = "boundary mismatch %.0f m" % st["boundary_mismatch"]
            out.append("%-4d %-8s %-10.2f %-10.2f %-10.0f %-10s %s" % (
                st["lev"], _g(st["dx"]), st["max_slope"], st["p99_slope"], st["relief"],
                "%.0f%%" % (100 * st["water_frac"]), extra))
    elif dem_name is None:
        out.append("")
        out.append("Terrain          : not checked (--dem none)")
    warns = list(plan.warnings) + list(twarn)
    if warns:
        out.append("")
        out.append("Warnings")
        for w in warns:
            out.append("  * " + w)
    return out


# --------------------------------------------------------------------------------------------

def namelist_wps(plan):
    f = plan.frame
    nd = len(plan.levels)
    ts = times(plan)
    start = wps_name(ts[0]) if ts else "YYYY-MM-DD_HH:00:00"
    end = wps_name(ts[-1]) if ts else "YYYY-MM-DD_HH:00:00"

    def row(vals):
        return ", ".join(str(v) for v in vals) + ","

    lines = [
        "! namelist.wps written by Exec/Tools/erfdomainadvisor.py for %s at %.4f, %.4f"
        % (plan.settings.phenomenon, plan.settings.lat, (plan.settings.lon + 180) % 360 - 180),
        "! d01 is ERF level 0; d02, d03, ... are ERF levels 1, 2, ... ERF reads each nest's place",
        "! from its met_em file (i_parent_start, j_parent_start, parent_grid_ratio).",
        "! Nests need only their first time: ERF reads one met_em file per refined level.",
        "&share",
        " wrf_core = 'ARW',",
        " max_dom = %d," % nd,
        " start_date = " + row(["'%s'" % start] * nd),
        " end_date   = " + row(["'%s'" % end] + ["'%s'" % start] * (nd - 1)),
        " interval_seconds = 3600,",
        " io_form_geogrid = 2,",
        "/",
        "",
        "&geogrid",
        " parent_id         = " + row([1] + list(range(1, nd))),
        " parent_grid_ratio = " + row([1] + [plan.ratio] * (nd - 1)),
        " i_parent_start    = " + row([1] + [lv.i_start + 1 for lv in plan.levels[1:]]),
        " j_parent_start    = " + row([1] + [lv.j_start + 1 for lv in plan.levels[1:]]),
        " e_we              = " + row([lv.nx + 1 for lv in plan.levels]),
        " e_sn              = " + row([lv.ny + 1 for lv in plan.levels]),
        " geog_data_res     = " + row(["'default'"] * nd),
        " dx = %s," % _exact(plan.base_dx),
        " dy = %s," % _exact(plan.base_dx),
        " map_proj  = 'lambert',",
        " ref_lat   = %.6f," % f.ref_lat,
        " ref_lon   = %.6f," % f.ref_lon,
        " truelat1  = %g," % f.truelat,
        " truelat2  = %g," % f.truelat,
        " stand_lon = %g," % f.stand_lon,
        " geog_data_path = 'USER_GEOG_DATA_PATH',  ! USER: the WPS static data directory",
        "/",
        "",
        "! ungrib: HRRR pressure-level GRIB2 files (wrfprsf) with the HRRR Vtable shipped with WPS",
        "! (ungrib/Variable_Tables/Vtable.raphrrr), linked as Vtable.",
        "&ungrib",
        " out_format = 'WPS',",
        " prefix = 'HRRR',",
        "/",
        "",
        "&metgrid",
        " fg_name = 'HRRR',",
        " io_form_metgrid = 2,",
        "/",
    ]
    return lines


# --------------------------------------------------------------------------------------------

def _closures(plan):
    pbl = ["None" if lv.regime == "LES" else "MYNN25" for lv in plan.levels]
    les = ["Smagorinsky" if lv.regime == "LES" else "None" for lv in plan.levels]
    return pbl, les


def _even_box(n, mult=1):
    """A max_grid_size that divides n into equal boxes: the divisor of n nearest 64 within
    [24, 128] that is a multiple of mult, or the multiple of mult nearest 64 when n has none.
    On a refined level mult is the ratio: AMReX chops there with max_grid_size / ratio."""
    divs = [d for d in range(24, 129) if n % d == 0 and d % mult == 0]
    return min(divs, key=lambda d: abs(d - 64)) if divs else mult * max(1, int(round(64.0 / mult)))


def inputs_deck(plan, command):
    s, f = plan.settings, plan.frame
    nlev = len(plan.levels) - 1
    lv0 = plan.levels[0]
    ts = times(plan)
    if ts:
        files0 = " ".join('"met_em.d01.%s.nc"' % wps_name(t) for t in ts)
        start_dt = ts[0].strftime("%Y-%m-%d %H:%M:%S")
    else:
        files0 = '"met_em.d01.FIRST_TIME.nc" "met_em.d01.SECOND_TIME.nc"'
        start_dt = None
    dt0 = lv0.dt
    steps_hour = int(round(3600.0 / dt0))
    pbl, les = _closures(plan)
    damp = damping_depth(plan.top)

    L = []
    a = L.append
    a("# ERF inputs written by Exec/Tools/erfdomainadvisor.py: a starting point, edit freely")
    a("# command: %s" % command)
    a("# phenomenon: %s (%s)" % (s.phenomenon, plan.title))
    a("# point: %.4f, %.4f at x = %.0f m, y = %.0f m" % (s.lat, (s.lon + 180) % 360 - 180,
                                                       plan.point_xy[0], plan.point_xy[1]))
    a("# Every value is a best guess; the report (report.txt) gives the reasons and warnings.")
    a("# Workflow: run WPS with namelist.wps on HRRR files, then ERF with this file in the")
    a("# directory of the met_em files. Check it first: python3 Exec/Tools/check_erf_inputs.py inputs")
    a("")
    a("erf.prob_name = \"Metgrid\"")
    a("")
    a("# ---- domain: must match namelist.wps d01 (nx = e_we - 1, ny = e_sn - 1, dx) ----")
    a("geometry.prob_lo     = 0.0 0.0 0.0")
    a("geometry.prob_hi     = %s %s %s" % (_exact(f.lx), _exact(f.ly), _exact(plan.top)))
    a("amr.n_cell           = %d %d %d   # dx = dy = %s m" % (lv0.nx, lv0.ny, plan.nz, _g(lv0.dx)))
    a("geometry.is_periodic = 0 0 0")
    if plan.stretch > 1.0:
        a("erf.initial_dz       = %s   # [m] first cell; ERF solves the stretching ratio (about %.4f)"
          % (_exact(plan.dz0), plan.stretch))
    else:
        a("# uniform vertical grid: dz = %s m" % _g(plan.top / plan.nz))
    a("")
    a("# ---- refinement: one WPS nest per level, placed by its met_em file ----")
    a("amr.max_level        = %d" % nlev)
    if nlev:
        a("amr.ref_ratio_vect   = %s   # x y z per level; z is not refined"
          % " ".join(["%d %d 1" % (plan.ratio, plan.ratio)] * nlev))
        for lv in plan.levels[1:]:
            a("# level %d: dx %s m, box %.1f-%.1f km by %.1f-%.1f km, %s"
              % (lv.lev, _g(lv.dx), lv.x0 / 1e3, lv.x1 / 1e3, lv.y0 / 1e3, lv.y1 / 1e3, lv.regime))
        a("erf.coupling_type    = OneWay   # coarse to fine only")
        a("#erf.coupling_type   = TwoWay   # fine levels feed back to their parents")
        a("erf.regrid_int       = -1       # the nests stay where WPS put them")
    a("")
    a("# ---- terrain ----")
    a("erf.terrain_type      = StaticFittedMesh   # required for Metgrid; terrain from the met_em files")
    a("erf.terrain_smoothing = 0   # basic terrain-following; ERF's Metgrid path allows only 0 when")
    a("                            # amr.max_level > 0")
    a("")
    a("# ---- boundaries ----")
    a("zlo.type = \"surface_layer\"")
    a("zhi.type = \"SlipWall\"")
    a("erf.real_width     = %d   # cells of the lateral relaxation zone" % s.real_width)
    a("")
    a("# ---- time ----")
    a("stop_time          = %s   # [s] %s h" % (_g(plan.hours * 3600), _g(plan.hours)))
    if start_dt:
        a("#start_datetime    = \"%s\"   # optional; must equal the first file's time" % start_dt)
    a("erf.fixed_dt       = %s   # [s] level 0; finer levels take %d steps per parent step"
      % (_g(dt0), plan.ratio))
    a("# erf.fixed_dt keeps the step-counted output below on the hour; with erf.cfl instead,")
    a("# replace the *_int keys by erf.plot_per_1, erf.plot2d_per_1 and erf.check_per [s].")
    a("#erf.cfl           = 0.8")
    a("erf.use_gravity    = true")
    a("")
    a("# ---- advection ----")
    a("erf.dycore_horiz_adv_type  = Upwind_5th")
    a("erf.dycore_vert_adv_type   = Upwind_3rd")
    a("erf.dryscal_horiz_adv_type = WENO5")
    a("erf.dryscal_vert_adv_type  = WENO5")
    a("")
    a("# ---- turbulence: per level; PBL scheme where dx >= 0.2 zi (zi = %g m), LES below ----" % plan.zi)
    a("erf.pbl_type        = %s" % " ".join(pbl))
    a("erf.les_type        = %s" % " ".join(les))
    if "Smagorinsky" in les:
        a("erf.Cs              = 0.25")
    a("erf.molec_diff_type = None")
    a("")
    a("# ---- surface layer ----")
    a("erf.surface_layer.flux_type = moeng")
    a("erf.most.z0 = 0.1   # [m] USER: about 2e-4 over water, 0.5 to 2 over cities; "
      "keep z0 below dz0 / 2 = %s m (raise --dz0 for a rougher surface)" % _g(plan.dz0 / 2))
    a("")
    a("# ---- Coriolis and damping ----")
    a("erf.use_coriolis   = true")
    a("erf.latitude       = %.4f" % s.lat)
    if max(f.lx, f.ly) > 500e3:
        a("erf.variable_coriolis = true   # Coriolis parameter from the file's latitudes")
    a("erf.rayleigh_damp_W  = true")
    a("erf.rayleigh_dampcoef = %s   # [1/s]; rate x dt kept at or below 0.25"
      % _g(min(0.2, 0.25 / dt0)))
    a("erf.rayleigh_zdamp   = %s   # [m] depth below the top" % _g(damp))
    a("")
    a("# ---- moisture ----")
    a("erf.moisture_model = %s" % plan.moisture)
    a("")
    a("# ---- initial and boundary data ----")
    a("erf.init_type      = Metgrid")
    a("erf.nc_init_file_0 = %s%s" % (files0, "" if ts else "   # USER: hourly met_em.d01 files"))
    for lv in plan.levels[1:]:
        name = ("met_em.d%02d.%s.nc" % (lv.lev + 1, wps_name(ts[0]))) if ts else \
            "met_em.d%02d.FIRST_TIME.nc" % (lv.lev + 1)
        a("erf.nc_init_file_%d = \"%s\"%s" % (lv.lev, name,
                                              "" if ts else "   # USER: first time of this nest"))
    a("")
    if plan.fire:
        fl = plan.levels[-1]
        a("# ---- fire: needs a build with a fire engine (in ERF-Hazard, FLARE) ----")
        a("erf.fire.enable       = true")
        a("#erf.fire.engine      = USER   # when more than one engine is built")
        a("erf.fire.anchor_level = -1     # the finest level, %d (dx %s m)" % (fl.lev, _g(fl.dx)))
        a("erf.fire.grid_ratio   = %d      # fire cells of %s m (--fire-dx %s m)"
          % (plan.fire_ratio, _g(fl.dx / plan.fire_ratio), _g(s.fire_dx)))
        a("# Every x,y box length on level %d must be a multiple of erf.fire.grid_ratio: its"
          % fl.lev)
        a("# %d x %d cells are multiples of %d, and max_grid_size %d cuts them into equal boxes."
          % (fl.nx, fl.ny, plan.fire_box, plan.fire_box))
        a("# Run on at most %d ranks, or AMReX cuts the boxes further." % (
            (fl.nx // plan.fire_box) * (fl.ny // plan.fire_box)))
        for d, attr in (("x", "nx"), ("y", "ny")):
            mgs = [_even_box(getattr(lv, attr), 1 if lv.lev == 0 else plan.ratio)
                   for lv in plan.levels[:-1]] + [plan.fire_box]
            a("amr.max_grid_size_%s   = %s" % (d, " ".join(map(str, mgs))))
        a("#erf.fire.smoke_enable = true   # smoke tracer; needs erf.transport_scalar = true")
        a("# USER: ignition, fuel map and fuel moisture are the engine's own keys; see its docs.")
        a("")
    a("# ---- output ----")
    a("erf.v            = 1")
    a("amr.v            = 1")
    a("erf.sum_interval = %d" % steps_hour)
    a("erf.data_log     = surf.log")
    a("erf.plot_file_1  = plt")
    a("erf.plot_int_1   = %d   # one per hour" % steps_hour)
    a("erf.plot_vars_1  = density x_velocity y_velocity z_velocity pressure theta temp rhoKE "
      "Kmv Khv z_phys" + (" qv qc" if plan.moisture != "None" else ""))
    a("erf.plot2d_file_1 = plt2d")
    a("erf.plot2d_int_1  = %d   # every 15 minutes" % max(1, steps_hour // 4))
    a("erf.plot2d_vars_1 = z_surf surf_pres u_star t_star Olen u_10m v_10m")
    a("erf.check_file   = chk")
    a("erf.check_int    = %d   # every 3 hours" % (3 * steps_hour))
    return L


def plan_json(plan, tstats=None, twarn=()):
    s, f = plan.settings, plan.frame
    d = {
        "settings": {k: getattr(s, k) for k in s.__dataclass_fields__},
        "title": plan.title,
        "target_dx": plan.target_dx, "base_dx": plan.base_dx, "ratio": plan.ratio,
        "domain": {"lx": f.lx, "ly": f.ly, "ref_lat": f.ref_lat, "ref_lon": f.ref_lon,
                   "truelat": f.truelat, "stand_lon": f.stand_lon},
        "point_xy": list(plan.point_xy),
        "vertical": {"top": plan.top, "nz": plan.nz, "dz0": plan.dz0, "stretch": plan.stretch,
                     "dz_top": plan.dz_top},
        "levels": [{k: getattr(lv, k) for k in lv.__dataclass_fields__} for lv in plan.levels],
        "cells": plan.cells, "core_hours": plan.core_hours, "memory_gb": plan.memory_gb,
        "hrrr_margin_cells": plan.hrrr_margin,
        "fire": {"enable": plan.fire, "grid_ratio": plan.fire_ratio, "box": plan.fire_box},
        "warnings": list(plan.warnings) + list(twarn),
        "terrain": tstats,
    }
    return json.dumps(d, indent=2)
