"""Command-line options and the per-setup defaults.

Every option left unset takes the default of the chosen ``--setup``.  The
resolved options are one ``argparse.Namespace`` (the "spec") that the
constraints, the recipes and the dummy-file writers read.
"""

import argparse
import json
import math

from .presets import PRESETS, preset_lines
from typing import List

CASES = ["ideal", "real"]
SETUPS = ["abl", "hill", "bubble", "density_current", "wrfinput", "metgrid"]
REAL_SETUPS = ["wrfinput", "metgrid"]
INITS = ["sounding", "uniform", "isentropic"]
DYCORES = ["compressible", "anelastic"]
TURBULENCE = ["none", "dns", "smagorinsky", "smagorinsky2d", "deardorff", "keqn",
              "myj", "mynn25", "mynnedmf", "ysu", "ysunew", "mrf"]
PBL = ["myj", "mynn25", "mynnedmf", "ysu", "ysunew", "mrf"]
TKE_CLOSURES = ["deardorff", "keqn", "myj", "mynn25", "mynnedmf"]
SURFACES = ["slip", "noslip", "most"]
ROUGHNESS = ["const", "file"]
SURFACE_HEAT = ["adiabatic", "flux", "temp", "file"]
TERRAIN = ["none", "analytic", "file", "file_usgs"]
HILLS = ["woa", "schar", "cos4"]
TERRAIN_METHODS = ["fitted", "immersed", "eb"]
BUILDINGS = ["none", "file"]
LATERAL = ["periodic", "slipwall", "outflow", "inflow_outflow", "open", "real"]
# the value of erf.moisture_model for each choice (SuperDroplets and SBM need their own builds and inputs)
MOISTURE = {"none": "None", "kessler": "Kessler", "kessler_norain": "Kessler_NoRain", "satadj": "SatAdj",
            "moistnocondensation": "MoistNoCondensation", "sam": "SAM", "sam_noice": "SAM_NoIce",
            "sam_noprecip_noice": "SAM_NoPrecip_NoIce", "morrison": "Morrison", "morrison_noice": "Morrison_NoIce",
            "wsm6": "WSM6", "wdm6": "WDM6"}
RADIATION = ["none", "twostream", "rrtmgp"]
LSM = ["none", "slm", "noahmp"]
CANOPY = ["none", "file"]
WINDFARM = ["none", "fitch", "ewp", "simplead"]
# refinement by field tests (ERF_Tagging.cpp): the fields a condition may name, and the tests
REFINE = ["box", "field"]
REFINE_FIELDS = ["density", "theta", "scalar", "qv", "qc", "qi", "qr", "qs", "qg", "qt", "vorticity", "helicity",
                 "max_reflectivity", "u", "v", "w", "omega", "terrain_blanking"]
REFINE_TESTS = {"greater": "value_greater", "less": "value_less", "gradient": "adjacent_difference_greater"}
REFINE_Z = ["above_ground", "physical"]
FORCING = ["none", "pressure_gradient", "geostrophic", "geo_table"]
DAMPING = ["none", "rayleigh", "sponge"]
OUTPUTS = ["plt", "plt2d", "netcdf", "chk", "logs", "profiles", "line", "plane", "stations", "subvol"]

# name, default domain (lo, hi), cells, and the options a setup starts from
SETUP_DEFAULTS = {
    "abl": dict(
        prob_name="ABL", lo=(0.0, 0.0, 0.0), hi=(1000.0, 1000.0, 1000.0), n=(64, 64, 64),
        init="sounding", turbulence="smagorinsky", surface="most", lateral="periodic",
        forcing="none", damping="none", wind=(10.0, 0.0), umax=10.0, stop_time=3600.0,
        terrain="none", about="atmospheric boundary layer; flat, or over terrain with --terrain"),
    "hill": dict(
        prob_name="Flow over Witch of Agnesi Hill", lo=(0.0, 0.0, 0.0), hi=(4000.0, 400.0, 2000.0), n=(80, 8, 40),
        init="isentropic", turbulence="smagorinsky", surface="noslip", lateral="periodic",
        forcing="none", damping="rayleigh", wind=(10.0, 0.0), umax=10.0, stop_time=1800.0,
        terrain="analytic", about="uniform flow over an analytic hill or ridge"),
    "bubble": dict(
        prob_name="Bubble", lo=(0.0, 0.0, 0.0), hi=(20000.0, 20000.0, 10000.0), n=(100, 100, 50),
        init="isentropic", turbulence="none", surface="slip", lateral="slipwall",
        forcing="none", damping="none", wind=(0.0, 0.0), umax=20.0, stop_time=1000.0,
        terrain="none", about="warm or cold bubble in a neutral atmosphere at rest"),
    "density_current": dict(
        prob_name="Density Current", lo=(-25600.0, 0.0, 0.0), hi=(25600.0, 100.0, 6400.0), n=(512, 1, 64),
        init="isentropic", turbulence="dns", surface="slip", lateral="outflow",
        forcing="none", damping="none", wind=(0.0, 0.0), umax=40.0, stop_time=900.0,
        terrain="none", about="cold bubble falling onto a flat surface (Straka et al. 1993 layout)"),
    "wrfinput": dict(
        prob_name="WPS", lo=(0.0, 0.0, 0.0), hi=(300000.0, 300000.0, 16000.0), n=(100, 100, 50),
        init="wrfinput", turbulence="mynn25", surface="most", lateral="real",
        forcing="none", damping="rayleigh", wind=(10.0, 0.0), umax=80.0, stop_time=3600.0,
        terrain="none", moisture="kessler", about="real-data case from WRF wrfinput and wrfbdy files"),
    "metgrid": dict(
        prob_name="Metgrid", lo=(0.0, 0.0, 0.0), hi=(300000.0, 300000.0, 16000.0), n=(100, 100, 50),
        init="metgrid", turbulence="mynn25", surface="most", lateral="real",
        forcing="none", damping="rayleigh", wind=(10.0, 0.0), umax=80.0, stop_time=3600.0,
        terrain="none", moisture="kessler", about="real-data case from WPS met_em files"),
}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="erfcreateinput.py",
        description="Write a commented ERF inputs file for a chosen configuration, with the "
                    "alternatives of each choice as comments, and the dummy auxiliary files the "
                    "deck names.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = ap.add_argument_group("presets and case files")
    g.add_argument("--preset", choices=sorted(PRESETS), default=None,
                   help="a named set of options (see --list); options given on the command line replace its values")
    g.add_argument("--case-file", default=None,
                   help="a JSON file of options, as written to case.json in every case directory; command-line "
                        "options replace its values")

    g = ap.add_argument_group("case")
    g.add_argument("--case", choices=CASES, default=None,
                   help="idealised or real-data case; follows --setup when unset")
    g.add_argument("--setup", choices=SETUPS, default="abl", help="problem setup (erf.prob_name)")
    g.add_argument("--init", choices=INITS, default=None, help="base state and initial profiles")
    g.add_argument("--dycore", choices=DYCORES, default="compressible", help="equation set")

    g = ap.add_argument_group("terrain, buildings and vertical grid")
    g.add_argument("--terrain", choices=TERRAIN, default=None,
                   help="terrain height: none, an analytic hill, an x-y-height table, or a USGS-style table")
    g.add_argument("--terrain-method", choices=TERRAIN_METHODS, default="fitted",
                   help="how terrain enters the grid: terrain-fitted mesh, immersed forcing, or embedded boundary")
    g.add_argument("--hill", choices=HILLS, default="woa",
                   help="shape for --terrain analytic: Witch of Agnesi, Schar mountain, cos^4 hill")
    g.add_argument("--hill-height", type=float, default=None, help="hill height [m]; default 5 percent of the domain depth")
    g.add_argument("--hill-width", type=float, default=None,
                   help="hill half-width [m]; default 5 percent of the x extent (25 percent for schar)")
    g.add_argument("--buildings", choices=BUILDINGS, default="none",
                   help="buildings by immersed forcing, from an x-y-height table")
    g.add_argument("--dz0", type=float, default=None,
                   help="first-cell height [m] of a vertically stretched grid; uniform cells when unset")

    g = ap.add_argument_group("physics")
    g.add_argument("--turbulence", choices=TURBULENCE, default=None, help="turbulence closure")
    g.add_argument("--surface", choices=SURFACES, default=None, help="lower boundary")
    g.add_argument("--roughness", choices=ROUGHNESS, default="const",
                   help="roughness length: one value, or a file of x y z0 (surface most)")
    g.add_argument("--surface-heat", choices=SURFACE_HEAT, default="adiabatic",
                   help="surface heating (surface most): none, a kinematic flux, a surface "
                        "temperature, or a time table file")
    g.add_argument("--heat-flux", type=float, default=0.1,
                   help="kinematic surface heat flux [K m/s] for --surface-heat flux; negative cools")
    g.add_argument("--moisture", choices=sorted(MOISTURE), default=None, help="moisture model")
    g.add_argument("--radiation", choices=RADIATION, default="none", help="radiation model")
    g.add_argument("--lsm", choices=LSM, default="none", help="land-surface model")
    g.add_argument("--canopy", choices=CANOPY, default="none", help="forest canopy drag from a table of tree patches")
    g.add_argument("--windfarm", choices=WINDFARM, default="none",
                   help="wind-farm parameterisation; needs a build with ERF_ENABLE_WINDFARM")
    g.add_argument("--amr", type=int, default=0, help="number of refined levels (0 to 2)")
    g.add_argument("--refine", choices=REFINE, default="box",
                   help="with --amr: static nested boxes at the domain centre, or cells tagged by field tests")
    g.add_argument("--refine-on", nargs="+", default=None, metavar="FIELD:TEST:VALUE",
                   help="with --refine field: one or more conditions, each a field (%s), a test (%s) and a value, "
                        "for example theta:less:299 w:greater:2" % (", ".join(REFINE_FIELDS), ", ".join(REFINE_TESTS)))
    g.add_argument("--refine-combine", choices=["and", "or"], default="and",
                   help="with several conditions: tag where all hold, or where any holds")
    g.add_argument("--refine-box", type=float, nargs="+", default=None, metavar="M",
                   help="with --refine field: restrict the tests to a region, as XLO YLO XHI YHI or "
                        "XLO YLO ZLO XHI YHI ZHI [m]")
    g.add_argument("--refine-z", choices=REFINE_Z, default=None,
                   help="how the z extent of --refine-box is measured; default above_ground on a fitted or flat "
                        "grid, physical with immersed or EB terrain")
    g.add_argument("--lateral", choices=LATERAL, default=None, help="lateral boundaries")
    g.add_argument("--forcing", choices=FORCING, default=None, help="large-scale driver")
    g.add_argument("--coriolis", action="store_true", help="turn the Coriolis force on")
    g.add_argument("--latitude", type=float, default=45.0, help="latitude [deg] for the Coriolis force")
    g.add_argument("--damping", choices=DAMPING, default=None, help="damping layer below the top")

    g = ap.add_argument_group("grid and time")
    g.add_argument("--nx", type=int, default=None, help="cells in x")
    g.add_argument("--ny", type=int, default=None, help="cells in y")
    g.add_argument("--nz", type=int, default=None, help="cells in z")
    g.add_argument("--domain", type=float, nargs=3, default=None, metavar=("LX", "LY", "LZ"),
                   help="domain extent [m]; the lower corner stays at the setup's default")
    g.add_argument("--dt", type=float, default=None, help="fixed time step [s]; adaptive when unset")
    g.add_argument("--stop-time", type=float, default=None, help="run length [s]")
    g.add_argument("--wind", type=float, nargs=2, default=None, metavar=("U", "V"),
                   help="initial and driving wind [m/s]")
    g.add_argument("--theta", type=float, default=300.0, help="reference potential temperature [K]")
    g.add_argument("--z0", type=float, default=0.1, help="roughness length [m]")

    g = ap.add_argument_group("output")
    g.add_argument("--output", default="plt,chk,logs",
                   help="comma-separated list of active outputs (%s); the others are written as "
                        "comments" % ", ".join(OUTPUTS))
    g.add_argument("--outdir", default="erf_case", help="directory to write the case into")
    g.add_argument("--name", default="inputs", help="name of the inputs file")
    g.add_argument("--force", action="store_true", help="overwrite files that already exist")
    g.add_argument("--list", action="store_true", help="print every choice and exit")
    g.add_argument("--self-test", action="store_true", help="run the built-in tests and exit")
    g.add_argument("--schema", default=None,
                   help="with --self-test: an inputs schema from build_erf_schema.py; every key the "
                        "generator can write must be in it")
    return ap


# options that describe the case; the others say where to write it or what to do
NOT_IN_CASE = ("outdir", "name", "force", "list", "self_test", "schema", "preset", "case_file")


class OptionError(Exception):
    """A preset or case file that cannot be applied."""


def parse(argv: List[str]) -> argparse.Namespace:
    """Parse the command line, then fill options not given there from --case-file and --preset."""
    parser = build_parser()
    args = parser.parse_args(argv)
    given = set()
    for action in parser._actions:
        for flag in action.option_strings:
            if any(a == flag or a.startswith(flag + "=") for a in argv):
                given.add(action.dest)
    known = {a.dest: a for a in parser._actions if a.dest not in NOT_IN_CASE and a.dest != "help"}
    layers = []
    if args.preset:
        layers.append(("--preset " + args.preset, {k: v for k, v in PRESETS[args.preset].items() if k != "about"}))
    if args.case_file:
        try:
            with open(args.case_file, encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, ValueError) as exc:
            raise OptionError("--case-file %s cannot be read: %s" % (args.case_file, exc))
        if not isinstance(data, dict):
            raise OptionError("--case-file %s must hold one JSON object of option names and values" % args.case_file)
        layers.append(("--case-file " + args.case_file, data))
    for source, values in layers:            # the case file is applied after the preset, so it wins
        for name, value in values.items():
            if name not in known:
                raise OptionError("%s: '%s' is not an option of this tool" % (source, name))
            action = known[name]
            if action.choices is not None and value is not None and value not in action.choices:
                raise OptionError("%s: %s = %r is not one of %s" % (source, name, value, ", ".join(map(str, action.choices))))
            if name not in given:
                setattr(args, name, value)
    return args


def case_options(args: argparse.Namespace) -> dict:
    """The options that describe the case, for case.json."""
    return {k: v for k, v in sorted(vars(args).items()) if k not in NOT_IN_CASE}


def resolve(args: argparse.Namespace) -> argparse.Namespace:
    """Fill unset options from the setup's defaults and derive the grid."""
    d = SETUP_DEFAULTS[args.setup]
    s = argparse.Namespace(**vars(args))
    s.real = s.setup in REAL_SETUPS
    s.case_given = s.case
    if s.case is None:
        s.case = "real" if s.real else "ideal"
    s.init_given = s.init
    if s.moisture is None:
        s.moisture = d.get("moisture", "none")
    if s.real:
        s.coriolis = True          # a real-data case is on the rotating Earth
    s.prob_name = d["prob_name"]
    for name in ("init", "turbulence", "surface", "lateral", "forcing", "damping", "wind", "stop_time", "terrain"):
        if getattr(s, name, None) is None:
            setattr(s, name, d[name])
    s.wind = tuple(s.wind)
    # largest speed expected in the run, for the time-step estimate; the buoyant setups start at rest
    s.umax = max(abs(s.wind[0]), abs(s.wind[1]), d["umax"])
    s.at_rest = s.wind == (0.0, 0.0)
    # a real-data run starts and ends on the file's clock; ideal runs on elapsed seconds
    s.lo = tuple(d["lo"])
    if s.domain is not None:
        s.hi = tuple(lo + ext for lo, ext in zip(s.lo, s.domain))
    else:
        s.hi = tuple(d["hi"])
    s.n = (s.nx if s.nx is not None else d["n"][0],
           s.ny if s.ny is not None else d["n"][1],
           s.nz if s.nz is not None else d["n"][2])
    s.dx = tuple((hi - lo) / n if n > 0 else 0.0 for lo, hi, n in zip(s.lo, s.hi, s.n))
    s.outputs = [o.strip() for o in s.output.split(",") if o.strip()]
    # periodic: x and y; slipwall: four walls; the others act on the x faces with y periodic
    s.periodic = {"periodic": (1, 1, 0), "slipwall": (0, 0, 0), "real": (0, 0, 0)}.get(s.lateral, (0, 1, 0))
    s.method = s.terrain_method if s.terrain != "none" else "none"
    if s.real:
        s.method = "fitted"        # the terrain comes from the NetCDF file; ERF requires the fitted mesh
    s.moist = s.moisture != "none"
    s.stretched = s.dz0 is not None
    s.buildings_on = s.buildings != "none"
    # with immersed terrain the wall law sits on the immersed surface, not on the zlo face
    s.zlo_most = s.surface == "most" and s.method != "immersed"
    lx, ly, lz = (s.hi[k] - s.lo[k] for k in range(3))
    if s.hill_height is None:
        s.hill_height = 0.05 * lz
    if s.hill_width is None:
        # the Schar ripples need the wider envelope to be resolved; the cos^4 hill ends at 4 L
        s.hill_width = {"schar": 0.25 * lx, "cos4": 0.05 * min(lx, ly) if s.n[1] >= RIDGE_BELOW else 0.05 * lx}.get(
            s.hill, 0.05 * lx)
    # with few cells across y the terrain is a ridge across x, uniform in y
    s.ridge = s.n[1] < RIDGE_BELOW
    s.gravity = s.init != "uniform"
    s.levels = s.amr
    s.conditions, s.refine_errors = refine_conditions(s)
    if s.refine_z is None:
        s.refine_z = "physical" if s.method in ("immersed", "eb") else "above_ground"
    # the vertical grid: first-cell height, and the ratio ERF will solve from erf.initial_dz
    s.stretch = stretch_ratio(s.dz0, lz, s.n[2]) if s.stretched and s.n[2] > 0 and lz > 0 else None
    s.dz_first = s.dz0 if s.stretched else s.dx[2]
    # the time step the deck will take, for the explicit damping rates: the fixed step, or the advective
    # estimate (the compressible adaptive step is smaller, so this is the larger of the two)
    positive = [d for d in s.dx if d > 0.0]
    s.dt_guess = float("%.2g" % (0.5 * min(positive) / s.umax)) if positive and s.umax > 0 else 0.0
    s.dt_est = s.dt if s.dt is not None else s.dt_guess
    # an adaptive anelastic step follows the wind, cfl x dx / |u|, and can exceed the estimate above
    s.adaptive_anelastic = s.dycore == "anelastic" and s.dt is None and not s.at_rest
    if s.adaptive_anelastic:
        speeds = [abs(w) / d for w, d in zip(s.wind, s.dx[:2]) if d > 0.0 and w != 0.0]
        s.dt_est = max(s.dt_est, float("%.2g" % (0.5 / sum(speeds)))) if speeds else s.dt_est
    s.is_pbl = s.turbulence in PBL
    s.has_tke = s.turbulence in TKE_CLOSURES
    return s


def refine_conditions(s):
    """The conditions of --refine field as (field, test, value) and the messages for those that cannot be read."""
    if s.amr <= 0 or s.refine != "field":
        return [], []
    given = s.refine_on
    if given is None:
        # a default that tags something in each setup; the user is asked to review it
        given = {"density_current": ["theta:less:%r" % (s.theta - 1.0)],
                 "bubble": ["theta:greater:%r" % (s.theta + 0.5)]}.get(s.setup, ["vorticity:greater:0.05"])
    conditions, errors = [], []
    for text in given:
        parts = str(text).split(":")
        if len(parts) != 3:
            errors.append("--refine-on %s: write each condition as FIELD:TEST:VALUE, for example theta:less:299" % text)
            continue
        field, test, value = parts
        try:
            number = float(value)
        except ValueError:
            number = float("nan")
        if field.lower() not in REFINE_FIELDS:
            errors.append("--refine-on %s: '%s' is not one of %s" % (text, field, ", ".join(REFINE_FIELDS)))
        elif test not in REFINE_TESTS:
            errors.append("--refine-on %s: '%s' is not one of %s" % (text, test, ", ".join(REFINE_TESTS)))
        elif not math.isfinite(number):
            errors.append("--refine-on %s: '%s' is not a finite number" % (text, value))
        else:
            conditions.append((field.lower(), test, number))
    return conditions, errors


RIDGE_BELOW = 16     # fewer cells than this across y: hills are written as ridges across x


def stretch_ratio(dz0: float, depth: float, nz: int):
    """The ratio ERF solves from erf.initial_dz alone, or None when its Newton iteration does not converge.

    The same iteration as geometric_stretch in Source/Utils/ERF_GridUtils.H: start at 1.03, at most 50
    steps, converged when the column depth is met to 1e-8 m.
    """
    tol = 1.0e-8
    if not (dz0 > 0.0 and depth > 0.0 and nz >= 1) or not all(math.isfinite(v) for v in (dz0, depth)):
        return None
    if nz == 1 or dz0 >= depth / nz:
        return 1.0

    def miss(r):
        return dz0 * (r ** nz - 1.0) / (r - 1.0) - depth

    try:
        r = 1.03
        f = miss(r)
        for _ in range(50):
            if abs(f) <= tol:
                break
            slope = dz0 * (nz * r ** (nz - 1) * (r - 1.0) - r ** nz + 1.0) / (r - 1.0) ** 2
            r = max(r - f / slope, 1.0 + tol)
            f = miss(r)
    except (OverflowError, ZeroDivisionError):
        return None
    return r if abs(f) <= tol else None


def list_choices() -> List[str]:
    out = ["--preset      :"] + preset_lines() + [""]
    out += ["--case        : " + ", ".join(CASES) + "   (real: --setup " + ", ".join(REAL_SETUPS) + ")",
            "--setup       :"]
    for name in SETUPS:
        d = SETUP_DEFAULTS[name]
        out.append("    %-16s %s" % (name, d["about"]))
        out.append("    %-16s defaults: init %s, turbulence %s, surface %s, lateral %s"
                   % ("", d["init"], d["turbulence"], d["surface"], d["lateral"]))
    out += ["--moisture    : " + ", ".join(sorted(MOISTURE)),
            "--radiation   : " + ", ".join(RADIATION),
            "--lsm         : " + ", ".join(LSM),
            "--canopy      : " + ", ".join(CANOPY),
            "--windfarm    : " + ", ".join(WINDFARM),
            "--amr         : 0, 1, 2",
            "--refine      : " + ", ".join(REFINE) + "   (--refine-on FIELD:TEST:VALUE ..., tests: "
            + ", ".join(REFINE_TESTS) + ")"]
    out += ["--terrain     : " + ", ".join(TERRAIN) + "   (--hill " + ", ".join(HILLS) + ")",
            "--terrain-method: " + ", ".join(TERRAIN_METHODS),
            "--buildings   : " + ", ".join(BUILDINGS),
            "--init        : " + ", ".join(INITS),
            "--dycore      : " + ", ".join(DYCORES),
            "--turbulence  : " + ", ".join(TURBULENCE),
            "--surface     : " + ", ".join(SURFACES),
            "--roughness   : " + ", ".join(ROUGHNESS),
            "--surface-heat: " + ", ".join(SURFACE_HEAT),
            "--lateral     : " + ", ".join(LATERAL),
            "--forcing     : " + ", ".join(FORCING),
            "--damping     : " + ", ".join(DAMPING),
            "--output      : " + ", ".join(OUTPUTS) + " (comma-separated)"]
    return out
