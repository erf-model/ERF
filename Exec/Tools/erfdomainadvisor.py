#!/usr/bin/env python3
"""Suggest an ERF domain, base resolution and refinement levels for a HRRR-driven case.

    python3 Exec/Tools/erfdomainadvisor.py --list
    python3 Exec/Tools/erfdomainadvisor.py --lat 34.05 --lon -118.53 --phenomenon wildfire \\
        --wind-from 45 --start 2025-01-07_12 --outdir palisades
    python3 Exec/Tools/erfdomainadvisor.py --lat 40.0 --lon -105.3 --phenomenon downslope_windstorm \\
        --dem none --outdir boulder
    python3 Exec/Tools/erfdomainadvisor.py --self-test

From a point and a phenomenon it prints a plan (domain size, base cell size, refined levels,
vertical grid, time step, cost) with checks of the terrain slope at each level, the distance to
the HRRR grid edge and the gray zone, and writes report.txt, namelist.wps (d01 plus one nest per
refined level), a starting ERF inputs file and plan.json.  Every number is a best guess; see
erf_domain_advisor/README.md.

It complements the other two input tools: erfcreateinput.py writes a deck for a chosen
configuration, check_erf_inputs.py checks any deck before a run.  Run the checker on the inputs
file this tool writes.
"""

import argparse
import os
import shlex
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from erf_domain_advisor import phenomena, writers  # noqa: E402
from erf_domain_advisor.plan import PlanError, Settings, build  # noqa: E402
from erf_domain_advisor.terrain import TerrainError, check, make_dem  # noqa: E402


def parser():
    p = argparse.ArgumentParser(
        description="Suggest an ERF domain, base resolution and refinement levels for a point and a "
                    "phenomenon, for a case driven by HRRR through WPS metgrid.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--lat", type=float, help="latitude of the point of interest [deg], -80 to 80")
    p.add_argument("--lon", type=float, help="longitude of the point of interest [deg], east positive, -180 to 360")
    p.add_argument("--phenomenon", choices=phenomena.names(), help="what the run is for (see --list)")
    g = p.add_argument_group("overrides of the phenomenon's best guesses")
    g.add_argument("--target-dx", type=float, help="finest cell size [m], 1 to 3000")
    g.add_argument("--ratio", type=int, default=3, help="refinement ratio between levels, 2 or 3")
    g.add_argument("--max-base-dx", type=float, default=1000.0,
                   help="largest base (level-0) cell size [m], 300 to 3000")
    g.add_argument("--outer", type=float, nargs=2, metavar=("LX", "LY"),
                   help="smallest outer domain [km], 1 to 5000 each")
    g.add_argument("--finest", type=float, nargs=2, metavar=("LX", "LY"),
                   help="finest box around the point [km], 1 to 5000 each")
    g.add_argument("--fetch", type=float, help="upwind fetch added to the finest box [km], 0 to 200")
    g.add_argument("--top", type=float, help="model top [km], 1 to 40")
    g.add_argument("--dz0", type=float, help="first-cell height [m], 0.5 to 200")
    g.add_argument("--hours", type=float, help="run length [h], 0.25 to 240")
    g.add_argument("--wind-from", type=float,
                   help="direction the wind comes from [deg, 0 = north, 90 = east; 0 to below 360]; places the "
                        "fetch upwind")
    g.add_argument("--start", help="first time, a whole hour, YYYY-MM-DD_HH; names the met_em files")
    fg = g.add_mutually_exclusive_group()
    fg.add_argument("--fire", dest="fire", action="store_const", const=True, default=None,
                   help="write an erf.fire.* block with the fire grid on the finest level "
                        "(default: on for wildfire, off otherwise)")
    fg.add_argument("--no-fire", dest="fire", action="store_const", const=False,
                   help="no erf.fire.* block, also for wildfire")
    g.add_argument("--fire-dx", type=float, default=10.0,
                   help="largest fire cell [m], 0.5 to 100; sets erf.fire.grid_ratio")
    g = p.add_argument_group("rules")
    g.add_argument("--real-width", type=int, default=5, help="erf.real_width [cells], 1 to 50")
    g.add_argument("--nest-margin", type=int, default=12,
                   help="smallest distance of a refined box to its parent's edge [parent cells], 4 to 100")
    g.add_argument("--hrrr-buffer", type=int, default=5,
                   help="HRRR cells between the relaxation zone and the first refined box, and "
                        "between the outer domain and the HRRR grid edge, 0 to 100")
    g.add_argument("--max-stretch", type=float, default=1.08, help="largest vertical stretching ratio, above 1 and at most 1.3")
    g.add_argument("--cost-per-cell-step", type=float, default=2.0e-6,
                   help="core-seconds per cell and step, for the cost estimate (a guess), 1e-9 to 1e-3")
    g = p.add_argument_group("terrain")
    g.add_argument("--dem", default="terrarium",
                   help="elevation source: terrarium (AWS Terrain Tiles, downloaded and cached), "
                        "an ESRI ASCII grid file in degrees, or none")
    g.add_argument("--slope-caution", type=float, default=0.3, help="slope that gets a note, above 0 and below --slope-limit")
    g.add_argument("--slope-limit", type=float, default=0.6, help="slope that gets a warning, at most 5")
    g.add_argument("--mismatch", type=float, default=150.0,
                   help="boundary terrain mismatch that gets a warning [m], above 0 and at most 5000")
    p.add_argument("--outdir", help="write report.txt, namelist.wps, inputs and plan.json here")
    p.add_argument("--force", action="store_true", help="overwrite files in --outdir")
    p.add_argument("--plot", action="store_true",
                   help="also write domain.png and slopes.png in --outdir (needs matplotlib)")
    p.add_argument("--list", action="store_true", help="list the phenomena and their best guesses")
    p.add_argument("--self-test", action="store_true", help="run the built-in tests (offline)")
    return p


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    args = parser().parse_args(argv)
    if args.self_test:
        from erf_domain_advisor import selftest
        return selftest.main()
    if args.list:
        print("\n".join(phenomena.table_lines()))
        return 0
    missing = [n for n in ("lat", "lon", "phenomenon") if getattr(args, n) is None]
    if missing:
        print("erfdomainadvisor: needs --%s (or --list)" % ", --".join(missing), file=sys.stderr)
        return 2
    try:
        if not (0.0 < args.slope_caution < args.slope_limit <= 5.0):
            raise PlanError("--slope-caution and --slope-limit need 0 < caution < limit <= 5, got "
                            "%g and %g" % (args.slope_caution, args.slope_limit))
        if not (0.0 < args.mismatch <= 5000.0):
            raise PlanError("--mismatch must lie in (0, 5000] m, got %g" % args.mismatch)
        s = Settings(lat=args.lat, lon=args.lon, phenomenon=args.phenomenon,
                     target_dx=args.target_dx, ratio=args.ratio, max_base_dx=args.max_base_dx,
                     outer=tuple(args.outer) if args.outer else None,
                     finest=tuple(args.finest) if args.finest else None,
                     fetch=args.fetch, top=args.top, dz0=args.dz0, hours=args.hours,
                     wind_from=args.wind_from, real_width=args.real_width,
                     nest_margin=args.nest_margin, hrrr_buffer=args.hrrr_buffer,
                     max_stretch=args.max_stretch, cost_per_cell_step=args.cost_per_cell_step,
                     start=args.start, fire=args.fire, fire_dx=args.fire_dx)
        if args.start:
            writers.parse_start(args.start)
        if args.plot and not args.outdir:
            raise PlanError("--plot needs --outdir")
        if args.plot:
            try:
                import matplotlib  # noqa: F401
            except ImportError:
                raise PlanError("--plot needs matplotlib")
        names = ("report.txt", "namelist.wps", "inputs", "plan.json")
        if args.plot:
            names += ("domain.png", "slopes.png")
        if args.outdir and not args.force:
            stale = [n for n in names if os.path.exists(os.path.join(args.outdir, n))]
            if stale:
                raise PlanError("%s already in %s; use --force" % (", ".join(stale), args.outdir))
        plan = build(s)
        dem = make_dem(args.dem)
        tstats, twarn, terrain = None, [], {}
        if dem is not None:
            print("erfdomainadvisor: terrain from %s" % dem.name, file=sys.stderr)
            tstats, twarn = check(plan, dem, args.slope_caution, args.slope_limit, args.mismatch,
                                  keep=terrain)
    except (PlanError, TerrainError, ValueError, OSError) as exc:
        print("erfdomainadvisor: %s" % exc, file=sys.stderr)
        return 2
    rep = writers.report(plan, tstats, twarn, dem.name if dem else None)
    print("\n".join(rep))
    if args.outdir:
        names = ("report.txt", "namelist.wps", "inputs", "plan.json")
        cmd = "erfdomainadvisor.py " + " ".join(shlex.quote(a) for a in argv)
        content = {
            "report.txt": rep,
            "namelist.wps": writers.namelist_wps(plan),
            "inputs": writers.inputs_deck(plan, cmd),
            "plan.json": [writers.plan_json(plan, tstats, twarn)],
        }
        try:
            os.makedirs(args.outdir, exist_ok=True)
            for n in names:
                with open(os.path.join(args.outdir, n), "w", newline="\n") as f:
                    f.write("\n".join(content[n]) + "\n")
            if args.plot:
                from erf_domain_advisor import plot
                names = names + tuple(plot.write(plan, terrain, args.slope_caution, args.slope_limit,
                                                 args.outdir, dem.name if dem else None))
        except OSError as exc:
            print("erfdomainadvisor: cannot write to --outdir %s: %s" % (args.outdir, exc), file=sys.stderr)
            return 2
        print("\nWrote %s in %s" % (", ".join(names), args.outdir))
        print("Next: run WPS with namelist.wps on HRRR files, then check the deck:")
        print("  python3 Exec/Tools/check_erf_inputs.py %s" % os.path.join(args.outdir, "inputs"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
