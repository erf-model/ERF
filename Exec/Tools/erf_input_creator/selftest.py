"""Built-in tests of the generator.

``python3 Exec/Tools/erfcreateinput.py --self-test`` writes a matrix of cases
into a temporary directory and checks, for each one, that

* every key the deck holds, active or commented, is a key ERF reads (its name
  is a string literal under ``Source/``; with ``--schema`` the full key must be
  in the schema);
* no key is active twice (the writer refuses it);
* the inputs checker reports no ERROR and no WARNING for the deck and its
  dummy files;
* each dummy file has the layout its reader expects for the case's grid;
* the deck parses back to the chosen options.

It also checks that each refused combination and each out-of-range option is
refused with a message naming the option.

With the environment variable ``ERF_CREATEINPUT_EXE`` set to an ERF
executable, every case whose files are all generated is also run for two
steps; the run must finish and must not list any input as unused.
``ERF_CREATEINPUT_FILTER`` restricts the runs to cases whose options hold
that text, and ``ERF_CREATEINPUT_NETCDF`` includes the NetCDF plotfile case.
"""

import json
import os
import re
import subprocess
import sys
import tempfile
from typing import List, Optional

from .auxfiles import dummy_building, dummy_hill
from .deckwriter import DeckText
from .enums import source_root
from .generate import CaseError, build, write
from .options import (DAMPING, FORCING, LATERAL, MOISTURE, OUTPUTS, SETUPS, SURFACE_HEAT, TERRAIN_METHODS,
                      TURBULENCE, WINDFARM, OptionError, case_options, parse, resolve)
from .presets import PRESETS

# Options of the cases written. The first block varies one option at a time from each setup's
# defaults; the second holds combinations whose keys interact.
CASES: List[List[str]] = [["--setup", name] for name in SETUPS]
_EXTRA = {"none": ["--surface", "slip"], "ysu": ["--surface-heat", "flux", "--heat-flux", "-0.01"]}
CASES += [["--setup", "abl", "--turbulence", t] + _EXTRA.get(t, []) for t in TURBULENCE]
CASES += [["--setup", "abl", "--surface", "slip"],
          ["--setup", "abl", "--surface", "noslip"],
          ["--setup", "abl", "--roughness", "file", "--nx", "12", "--ny", "10"]]
CASES += [["--setup", "abl", "--surface-heat", h] for h in SURFACE_HEAT]
CASES += [["--setup", "abl", "--lateral", lat] + (["--wind", "0", "0"] if lat == "slipwall" else [])
          for lat in LATERAL if lat != "real"]
CASES += [["--setup", "abl", "--forcing", f] + (["--coriolis"] if f.startswith("geo") else []) for f in FORCING]
CASES += [["--setup", "abl", "--damping", d] for d in DAMPING]
CASES += [["--setup", "abl", "--init", "uniform"],
          ["--setup", "abl", "--init", "isentropic"],
          ["--setup", "abl", "--dycore", "anelastic"],
          ["--setup", "abl", "--dt", "0.5"],
          ["--setup", "abl", "--dycore", "anelastic", "--dt", "0.5"],
          ["--setup", "abl", "--output", ",".join(o for o in OUTPUTS if o != "netcdf")],
          ["--setup", "abl", "--output", "plt,netcdf"],
          ["--setup", "abl", "--output", ""],
          ["--setup", "abl", "--domain", "2000", "1000", "500", "--nx", "20", "--ny", "10", "--nz", "8",
           "--roughness", "file", "--surface-heat", "flux", "--forcing", "geo_table", "--coriolis",
           "--turbulence", "mynn25", "--damping", "rayleigh", "--output", "plt,plt2d,chk,logs,profiles,stations"],
          ["--setup", "abl", "--surface-heat", "file", "--theta", "290"],
          ["--setup", "bubble", "--dycore", "anelastic"],
          ["--setup", "bubble", "--init", "sounding"],
          ["--setup", "density_current", "--dycore", "anelastic"],
          ["--setup", "density_current", "--turbulence", "smagorinsky"]]

# terrain, immersed forcing, embedded boundary, buildings, vertical stretching
CASES += [["--setup", "hill", "--terrain-method", m] for m in TERRAIN_METHODS if m != "eb"]
CASES += [["--setup", "hill", "--hill", "woa"],
          ["--setup", "hill", "--hill", "schar", "--hill-width", "500"],
          ["--setup", "hill", "--hill", "cos4", "--hill-width", "40", "--domain", "1600", "400", "800",
           "--nx", "32", "--ny", "16", "--nz", "16"]]
CASES += [["--setup", "hill", "--dycore", "anelastic", "--terrain-method", m] for m in TERRAIN_METHODS]
CASES += [["--setup", "abl", "--terrain", "file", "--terrain-method", m] for m in TERRAIN_METHODS if m != "eb"]
CASES += [["--setup", "abl", "--terrain", "analytic", "--terrain-method", "eb", "--dycore", "anelastic",
           "--domain", "1600", "400", "800", "--nx", "32", "--ny", "8", "--nz", "16"],
          ["--setup", "hill", "--terrain", "file", "--terrain-method", "eb", "--dycore", "anelastic"]]
CASES += [["--setup", "abl", "--terrain", "file", "--buildings", "file", "--nx", "12", "--ny", "8", "--nz", "16"],
          ["--setup", "abl", "--terrain", "file_usgs"],
          ["--setup", "abl", "--terrain", "file_usgs", "--terrain-method", "immersed", "--nx", "12", "--ny", "10"],
          ["--setup", "abl", "--terrain", "analytic", "--hill", "cos4", "--hill-height", "80", "--hill-width", "60"],
          ["--setup", "abl", "--terrain", "file", "--roughness", "file", "--surface-heat", "flux"],
          ["--setup", "abl", "--terrain", "file", "--terrain-method", "immersed", "--surface-heat", "flux"],
          ["--setup", "abl", "--terrain", "file", "--terrain-method", "immersed", "--surface-heat", "temp"],
          ["--setup", "abl", "--terrain", "file", "--terrain-method", "immersed", "--surface", "noslip"],
          ["--setup", "hill", "--terrain-method", "eb", "--dycore", "anelastic", "--surface", "most", "--z0", "0.05"],
          ["--setup", "abl", "--dz0", "20"],
          ["--setup", "abl", "--dz0", "20", "--dycore", "anelastic"],
          ["--setup", "abl", "--terrain", "file", "--dz0", "20", "--turbulence", "mynn25"],
          ["--setup", "abl", "--buildings", "file"],
          ["--setup", "abl", "--buildings", "file", "--dycore", "anelastic"],
          ["--setup", "abl", "--terrain", "file", "--terrain-method", "immersed", "--buildings", "file"],
          ["--setup", "abl", "--terrain", "file", "--buildings", "file"],
          ["--setup", "hill", "--init", "sounding", "--ny", "1", "--nx", "32", "--nz", "16"]]

# moisture, radiation, land surface, canopy, wind farm, refinement, sponge
CASES += [["--setup", "abl", "--moisture", m] for m in sorted(MOISTURE)]
CASES += [["--setup", "abl", "--windfarm", w] for w in WINDFARM if w != "none"]
CASES += [["--setup", "abl", "--radiation", "twostream"],
          ["--setup", "abl", "--radiation", "twostream", "--moisture", "kessler", "--lsm", "slm"],
          ["--setup", "abl", "--radiation", "rrtmgp", "--moisture", "kessler"],
          ["--setup", "abl", "--moisture", "kessler", "--lsm", "slm"],
          ["--setup", "abl", "--moisture", "kessler", "--lsm", "noahmp"],
          ["--setup", "abl", "--canopy", "file"],
          ["--setup", "abl", "--canopy", "file", "--turbulence", "keqn", "--terrain", "file"],
          ["--setup", "abl", "--amr", "1"],
          ["--setup", "abl", "--amr", "2", "--dycore", "anelastic"],
          ["--setup", "abl", "--amr", "1", "--terrain", "file"],
          ["--setup", "abl", "--amr", "1", "--terrain", "file", "--terrain-method", "immersed"],
          ["--setup", "abl", "--damping", "sponge"],
          ["--setup", "hill", "--damping", "sponge", "--moisture", "kessler"],
          ["--setup", "hill", "--terrain-method", "eb", "--dycore", "anelastic", "--moisture", "kessler"],
          ["--setup", "bubble", "--moisture", "kessler"],
          # real-data cases
          ["--setup", "wrfinput"],
          ["--setup", "wrfinput", "--case", "real", "--moisture", "morrison", "--lsm", "noahmp", "--radiation", "rrtmgp",
           "--turbulence", "mrf", "--output", "plt,plt2d,chk,stations"],
          ["--setup", "metgrid"],
          ["--setup", "metgrid", "--moisture", "none", "--turbulence", "smagorinsky", "--dycore", "anelastic"]]

CASES += [["--preset", name] for name in sorted(PRESETS)]
CASES += [["--setup", "abl", "--theta", "290.12345", "--wind", "8", "3", "--z0", "0.03", "--surface-heat", "temp",
           "--latitude", "-33.5", "--coriolis", "--stop-time", "7200", "--name", "inputs_abl", "--heat-flux", "0.25"],
          ["--setup", "abl", "--init", "isentropic", "--theta", "290", "--wind", "8", "3"],
          ["--setup", "abl", "--init", "uniform", "--theta", "285", "--wind", "-6", "2"],
          ["--setup", "bubble", "--wind", "10", "0", "--lateral", "periodic"],
          ["--setup", "abl", "--surface-heat", "flux", "--heat-flux", "0.25"],
          ["--setup", "abl", "--forcing", "geo_table", "--coriolis", "--domain", "1000", "1000", "1000.0004"],
          ["--setup", "abl", "--domain", "1000", "1000", "333.3333333", "--nz", "12"],
          ["--setup", "abl", "--moisture", "sam", "--dz0", "20"],
          ["--setup", "abl", "--dycore", "anelastic", "--dz0", "20", "--damping", "sponge", "--dt", "3"],
          ["--setup", "abl", "--buildings", "file", "--canopy", "file"],
          ["--setup", "density_current", "--dycore", "anelastic", "--damping", "rayleigh"]]

# (options, text the message must hold)
REFUSED = [
    (["--surface", "most", "--turbulence", "none"], "--surface most"),
    (["--turbulence", "ysu", "--surface", "slip"], "--turbulence ysu"),
    (["--turbulence", "mynn25", "--surface", "noslip"], "--turbulence mynn25"),
    (["--turbulence", "ysu"], "--turbulence ysu"),
    (["--turbulence", "ysu", "--surface-heat", "flux"], "--heat-flux"),
    (["--forcing", "geostrophic"], "--coriolis"),
    (["--forcing", "geo_table"], "--coriolis"),
    (["--roughness", "file", "--surface", "slip"], "--roughness file"),
    (["--surface-heat", "flux", "--surface", "slip"], "--surface-heat flux"),
    (["--surface-heat", "temp", "--init", "uniform"], "--surface-heat temp"),
    (["--setup", "bubble", "--init", "uniform"], "--init uniform"),
    (["--setup", "density_current", "--init", "uniform"], "--init uniform"),
    (["--roughness", "file", "--nx", "1"], "--roughness file"),
    (["--terrain-method", "eb"], "needs terrain"),
    (["--setup", "hill", "--terrain-method", "eb"], "--terrain-method eb"),
    (["--terrain", "file", "--terrain-method", "eb", "--dycore", "anelastic", "--ny", "8"], "--setup hill only"),
    (["--terrain", "file", "--terrain-method", "eb", "--dycore", "anelastic", "--surface", "slip"], "--terrain-method eb"),
    (["--setup", "hill", "--terrain-method", "eb", "--dycore", "anelastic", "--canopy", "file"], "--canopy"),
    (["--setup", "hill", "--terrain-method", "eb", "--dycore", "anelastic", "--amr", "1"], "--amr"),
    (["--terrain", "file", "--terrain-method", "immersed", "--canopy", "file"], "--canopy"),
    (["--z0", "40"], "--z0"),
    (["--dz0", "0.15", "--nz", "16"], "--dz0"),
    (["--dz0", "0.2", "--z0", "0.1"], "--z0"),
    (["--turbulence", "keqn", "--surface", "slip"], "--turbulence keqn"),
    (["--turbulence", "none", "--surface", "slip", "--output", "plt,profiles"], "--output profiles"),
    (["--amr", "1", "--output", "plt,plane"], "--output plane"),
    (["--lateral", "slipwall"], "--lateral slipwall"),
    (["--lateral", "inflow_outflow", "--wind", "-5", "0"], "--lateral inflow_outflow"),
    (["--setup", "density_current", "--wind", "5", "0"], "--setup density_current"),
    (["--lsm", "slm", "--moisture", "kessler", "--surface", "slip", "--turbulence", "none"], "--lsm slm"),
    (["--surface-heat", "file", "--roughness", "file"], "--surface-heat file"),
    (["--setup", "hill", "--hill", "schar", "--hill-width", "100"], "--hill schar"),
    (["--setup", "hill", "--hill", "cos4"], "--hill cos4"),
    (["--terrain", "analytic", "--hill", "cos4", "--hill-width", "400"], "--hill cos4"),
    (["--heat-flux", "nan"], "--heat-flux"),
    (["--wind", "inf", "0"], "--wind"),
    (["--domain", "inf", "1", "1"], "--domain"),
    (["--nz", "1"], "--nz"),
    (["--name", "input_sounding"], "--name"),
    (["--name", "FILES_TO_REPLACE.txt"], "--name"),
    (["--name", "."], "--name"),
    (["--setup", "hill", "--hill-height", "1500"], "--hill-height"),
    (["--case", "real"], "--case real"),
    (["--setup", "wrfinput", "--case", "ideal"], "--case ideal"),
    (["--setup", "wrfinput", "--init", "sounding"], "--init"),
    (["--setup", "wrfinput", "--lateral", "periodic"], "--lateral real"),
    (["--lateral", "real"], "--lateral real"),
    (["--setup", "metgrid", "--terrain", "file"], "NetCDF file"),
    (["--setup", "metgrid", "--dz0", "50"], "NetCDF file"),
    (["--setup", "wrfinput", "--moisture", "none"], "--setup wrfinput"),
    (["--setup", "wrfinput", "--forcing", "pressure_gradient"], "--forcing"),
    (["--moisture", "sam", "--terrain", "file"], "--moisture sam"),
    (["--moisture", "sam_noice", "--terrain", "analytic"], "--moisture sam_noice"),
    (["--setup", "wrfinput", "--moisture", "sam"], "--moisture sam"),
    (["--lsm", "slm"], "--lsm slm"),
    (["--amr", "3"], "--amr"),
    (["--amr", "-1"], "--amr"),
    (["--amr", "1", "--nx", "18"], "--amr"),
    (["--setup", "wrfinput", "--damping", "sponge"], "--damping sponge"),
    (["--terrain", "file", "--terrain-method", "eb", "--turbulence", "deardorff", "--surface", "slip"], "--terrain-method eb"),
    (["--terrain", "file", "--terrain-method", "eb", "--roughness", "file"], "--roughness file"),
    (["--terrain", "file", "--terrain-method", "eb", "--buildings", "file"], "--buildings file"),
    (["--terrain", "file", "--terrain-method", "immersed", "--dz0", "5"], "--dz0"),
    (["--buildings", "file", "--dz0", "5"], "--dz0"),
    (["--terrain", "file", "--terrain-method", "immersed", "--roughness", "file"], "--roughness file"),
    (["--terrain", "file", "--terrain-method", "immersed", "--turbulence", "mrf"], "--turbulence mrf"),
    (["--terrain", "file_usgs", "--ny", "1"], "file_usgs"),
    (["--dz0", "0"], "--dz0"),
    (["--dz0", "500"], "--dz0"),
    (["--setup", "hill", "--hill-height", "-5"], "--hill-height"),
    (["--setup", "hill", "--hill-height", "5000"], "--hill-height"),
    (["--setup", "hill", "--hill-width", "0"], "--hill-width"),
    (["--output", "netcdf"], "--output netcdf"),
    (["--output", "plt,movie"], "'movie'"),
    (["--nx", "0"], "--nx"),
    (["--nz", "-4"], "--nz"),
    (["--domain", "100", "0", "100"], "--domain"),
    (["--dt", "0"], "--dt"),
    (["--stop-time", "-1"], "--stop-time"),
    (["--theta", "0"], "--theta"),
    (["--z0", "0"], "--z0"),
    (["--latitude", "91"], "--latitude"),
    (["--name", "a/b"], "--name"),
]


def _spec(options: List[str], outdir: str):
    return resolve(parse(options + ["--outdir", outdir]))


def _source_literals():
    """Every identifier-like string literal under Source/ and Exec/*.cpp: the pool of key names ERF reads."""
    root = source_root()
    names = set()
    literal = re.compile(r'"([A-Za-z_][A-Za-z0-9_.]*)"')
    tops = [os.path.join(root, "Source")]
    files = [os.path.join(root, "Exec", fn) for fn in ("ERF_Prob.cpp", "ERF_Prob.H")]
    for top in tops:
        for d, _, fns in os.walk(top):
            files += [os.path.join(d, fn) for fn in fns if fn.endswith((".H", ".cpp"))]
    for path in files:
        try:
            with open(path, errors="replace") as fh:
                for m in literal.finditer(fh.read()):
                    names.add(m.group(1).rsplit(".", 1)[-1])
        except OSError:
            continue
    return names


def _schema_knows(schema: dict, key: str) -> bool:
    if key in schema["keys"]:
        return True
    # keys whose prefix or suffix is built at run time (faces, stations, level numbers) are held as
    # patterns such as <face>.type or erf.nc_init_file_<n>
    for pattern in schema.get("patterns") or {}:
        regex = "".join(r"[A-Za-z0-9_]+" if part.startswith("<") else re.escape(part)
                        for part in re.split(r"(<[^>]*>)", str(pattern)))
        if re.fullmatch(regex, key) or re.fullmatch(regex, key.split(".", 1)[-1]):
            return True
    return False


_NAMED_PREFIX = re.compile(r"^(erf\.(mast\d+|box\d+)|[xyz](lo|hi))\.")


def _check_keys(deck, schema: Optional[dict], read_names, documented, failures, label) -> None:
    """Every key, active or commented, must be one ERF reads, under the prefix ERF reads it with."""
    by_leaf = {}
    for full in documented:
        by_leaf.setdefault(full.rsplit(".", 1)[-1], set()).add(full.rsplit(".", 1)[0] if "." in full else "")
    for key in deck.all_keys():
        leaf = key.rsplit(".", 1)[-1]
        prefix = key.rsplit(".", 1)[0] if "." in key else ""
        if key in ("max_step", "stop_time", "start_datetime", "stop_datetime"):
            continue
        # not string literals in ERF's Source/: the level number is appended at run time, and
        # ref_ratio_vect is read by AMReX's AmrMesh. Both are checked against the schema when one is given.
        composed = key in ("erf.nc_init_file_0", "amr.ref_ratio_vect")
        # the prefix: a key the docs list is right as written; a name the docs list only under other
        # prefixes is wrong here (erf.Cs written as amr.Cs). Keys under a user-chosen name (a station,
        # a refinement box, a face) and prob.* keys are not listed with their prefix in the docs.
        if key not in documented and not _NAMED_PREFIX.match(key) and not key.startswith("prob."):
            elsewhere = by_leaf.get(leaf, set())
            if elsewhere and prefix not in elsewhere:
                failures.append("%s: the key %s is not documented; the docs list %s under %s"
                                % (label, key, leaf, ", ".join(sorted(elsewhere))))
        # read by SLM through a helper that takes the name as an argument (get_layer_prop in ERF_SLM.cpp), which
        # the schema scanner does not follow; the names are string literals in Source/, and the ERF runs read them
        if key in ("slm.clay0", "slm.sand0", "slm.sw0", "slm.st0"):
            if leaf not in read_names:
                failures.append("%s: the key %s is not read anywhere under Source/" % (label, key))
            continue
        if read_names and leaf not in read_names and not composed:
            failures.append("%s: the key %s is not read anywhere under Source/" % (label, key))
        if schema is not None and not key.startswith("prob.") and not _schema_knows(schema, key):
            failures.append("%s: the key %s is not in the schema" % (label, key))


def _check_aux(s, outdir: str, failures, label) -> None:
    def rows(name):
        with open(os.path.join(outdir, name)) as fh:
            return [ln.split() for ln in fh.read().splitlines()]

    if s.init == "sounding":
        r = rows("input_sounding")
        if len(r[0]) != 3 or any(len(x) != 5 for x in r[1:]):
            failures.append("%s: input_sounding needs 3 values on line 1 and 5 on the others" % label)
        z = [float(x[0]) for x in r[1:]]
        if any(b <= a for a, b in zip(z, z[1:])) or z[0] != s.lo[2] or z[-1] < s.hi[2]:
            failures.append("%s: input_sounding heights must increase from the bottom to the domain top" % label)
        if r[1][1] != r[0][1] or r[1][2] != r[0][2]:
            failures.append("%s: input_sounding bottom row must repeat the surface theta and qv" % label)
    nxn, nyn = s.n[0] + 1, s.n[1] + 1
    for name, used, periodic in (("terrain.txt", s.terrain == "file", True), ("buildings.txt", s.buildings_on, False)):
        if not used:
            continue
        v = [float(x[0]) for x in rows(name) if len(x) == 1]
        if len(v) != 2 + nxn + nyn + nxn * nyn or (v[0], v[1]) != (nxn, nyn):
            failures.append("%s: %s needs nx, ny, %d x, %d y and %d heights, one per line" % (label, name, nxn, nyn, nxn * nyn))
            continue
        xs, ys, h = v[2:2 + nxn], v[2 + nxn:2 + nxn + nyn], v[2 + nxn + nyn:]
        if xs[0] != s.lo[0] or abs(xs[-1] - s.hi[0]) > 1e-6 or ys[0] != s.lo[1] or abs(ys[-1] - s.hi[1]) > 1e-6:
            failures.append("%s: %s does not span the domain" % (label, name))
        if max(h) <= 0.0 or min(h) < 0.0 or max(h) >= s.hi[2] - s.lo[2]:
            failures.append("%s: %s heights must lie in [0, domain depth) with a non-zero maximum" % (label, name))
        # y fastest: the first ny values belong to x[0], a lateral face, where the dummy terrain is zero
        if periodic and (any(h[j] != 0.0 for j in range(nyn)) or any(h[(nxn - 1) * nyn + j] != 0.0 for j in range(nyn))
                         or (nyn > 4 and any(h[i * nyn] != h[i * nyn + nyn - 1] for i in range(nxn)))):
            failures.append("%s: %s must be equal on opposite lateral faces, with y fastest" % (label, name))
        # y fastest: value n belongs to (x[n // ny], y[n % ny]); a grid with nx != ny tells the two orders apart
        shape = dummy_hill if name == "terrain.txt" else dummy_building
        if any(abs(h[i * nyn + j] - shape(s, xs[i], ys[j])) > 1e-5 for i in range(nxn) for j in range(nyn)):
            failures.append("%s: %s heights are not in y-fastest order" % (label, name))
    if s.terrain == "file_usgs":
        r = rows("terrain_usgs.txt")
        body = r[2:]
        if len(r[0]) != 2 or [int(x) for x in r[1]] != [nxn, nyn] or len(body) != nxn * nyn or any(len(x) != 3 for x in body):
            failures.append("%s: terrain_usgs.txt needs 'lon lat', 'nx ny' and %d rows of x y height" % (label, nxn * nyn))
        elif float(body[1][0]) <= float(body[0][0]) or body[1][1] != body[0][1] or float(body[nxn][1]) <= float(body[0][1]):
            failures.append("%s: terrain_usgs.txt must run x fastest" % label)
    if s.canopy == "file":
        r = rows("forest.txt")
        if not r or any(len(x) != 8 or x[0] not in ("1", "2") for x in r):
            failures.append("%s: forest.txt needs rows of 8 numbers starting with tree type 1 or 2" % label)
        elif not (s.lo[0] < float(r[0][1]) < s.hi[0] and s.lo[1] < float(r[0][2]) < s.hi[1]):
            failures.append("%s: forest.txt patch centre is outside the domain" % label)
    if s.windfarm != "none":
        r = rows("windturbines_loc_x_y.txt")
        if not r or any(len(x) != 2 for x in r) or not s.lo[0] < float(r[0][0]) < s.hi[0]:
            failures.append("%s: windturbines_loc_x_y.txt needs rows 'x y' inside the domain" % label)
        r = rows("windturbines_spec.tbl")
        ok = len(r) >= 4 and len(r[0]) == 1 and len(r[1]) == 4 and all(len(x) == 3 for x in r[2:])
        if not ok or int(r[0][0]) != len(r) - 2 or int(r[0][0]) < 2:
            failures.append("%s: windturbines_spec.tbl needs a row count, 4 turbine values and that many 3-column rows" % label)
        elif 0.5 * float(r[1][1]) > float(r[1][0]) or any(float(x[1]) > 1.0 for x in r[2:]):
            failures.append("%s: windturbines_spec.tbl needs rotor radius <= hub height and Ct <= 1" % label)
    if s.zlo_most and s.roughness == "file":
        r = rows("roughness.txt")
        nx, ny = s.n[0], s.n[1]
        if len(r) != (nx + 1) * (ny + 1) or any(len(x) != 3 for x in r):
            failures.append("%s: roughness.txt needs %d rows of x y z0" % (label, (nx + 1) * (ny + 1)))
        else:
            for j in (0, ny):
                for i in (0, 1, nx):
                    x, y, z0 = (float(v) for v in r[i + j * (nx + 1)])
                    if abs(x - (s.lo[0] + i * s.dx[0])) > 1e-4 or abs(y - (s.lo[1] + j * s.dx[1])) > 1e-4 or z0 <= 0:
                        failures.append("%s: roughness.txt row %d is not node (%d, %d)" % (label, i + j * (nx + 1), i, j))
    if s.zlo_most and s.surface_heat == "file":
        r = rows("sfc.txt")
        if len(r) < 3 or any(len(x) != 2 for x in r[1:]):
            failures.append("%s: sfc.txt needs a header and rows of day sst" % label)
    if s.forcing == "geo_table":
        r = [x for x in rows("geo_wind_table.txt") if x and not x[0].startswith("#")]
        z = [float(x[0]) for x in r]
        if any(len(x) != 3 for x in r) or z[0] > s.lo[2] or z[-1] < s.hi[2]:
            failures.append("%s: geo_wind_table.txt must span the domain with rows z Ug Vg" % label)


def _check_deck(s, outdir: str, failures, label, needed=()) -> None:
    """The checker's verdict, and the deck read back against the options."""
    from erf_input_checker import Report, RunContext, Severity, parse_deck, run_all_checks
    from erf_input_checker.deck import DeckError
    path = os.path.join(outdir, s.name)
    try:
        deck = parse_deck(path)
        deck.get_ints("amr.n_cell")
    except DeckError as exc:
        failures.append("%s: the deck does not parse: %s" % (label, exc))
        return
    report = Report(path)
    run_all_checks(deck, RunContext(deck), report)
    for f in report.sorted():
        # a file the tool cannot generate is reported missing; that is the expected state of a new case
        if f.key in needed and "not found" in f.message:
            continue
        if f.severity >= Severity.WARNING:
            failures.append("%s: checker %s" % (label, f.format().replace("\n", " ")))
    for key in deck.keys():
        if any("#" in word for word in deck.raw(key) or []):
            failures.append("%s: a comment is glued to the value of %s; ParmParse reads it as part of the value"
                            % (label, key))
    if (s.dycore == "anelastic") != (deck.get_str("erf.use_fft") == "true"):
        failures.append("%s: erf.use_fft = true must be active for anelastic cases, and only for them" % label)
    want = {"erf.compressibility_type": "FullyCompressible" if s.dycore == "compressible" else "Anelastic",
            "erf.use_gravity": "true" if s.gravity else "false",
            "zlo.type": "NoSlipWall" if s.method == "immersed" else
            {"most": "surface_layer", "slip": "SlipWall", "noslip": "NoSlipWall"}[s.surface],
            "erf.terrain_type": {"none": "None", "fitted": "StaticFittedMesh", "immersed": "ImmersedForcing",
                                 "eb": "EB"}[s.method],
            "erf.init_type": {"sounding": "Input_Sounding", "uniform": "Uniform", "isentropic": "Isentropic",
                              "wrfinput": "WRFInput", "metgrid": "Metgrid"}[s.init],
            "erf.moisture_model": MOISTURE[s.moisture]}
    for key, value in want.items():
        if deck.get_str(key) != value:
            failures.append("%s: %s is %s, expected %s" % (label, key, deck.get_str(key), value))
    def close(key, value):
        got = deck.get_real(key)
        if got is None or abs(got - value) > 1e-12 * max(1.0, abs(value)):
            failures.append("%s: %s is %s, expected %s" % (label, key, got, value))

    close("stop_time", s.stop_time)
    if s.init in ("uniform", "isentropic") and s.setup in ("abl", "hill"):
        close("prob.T_0", s.theta)
        close("prob.U_0", s.wind[0])
        close("prob.V_0", s.wind[1])
    if s.setup == "bubble" and not s.at_rest:
        close("prob.U_0", s.wind[0])
    if s.init == "sounding":
        with open(os.path.join(outdir, "input_sounding")) as fh:
            rows = [ln.split() for ln in fh.read().splitlines()]
        if float(rows[0][1]) != s.theta or [float(v) for v in rows[-1][3:5]] != list(s.wind):
            failures.append("%s: the sounding does not hold --theta and --wind" % label)
        if s.zlo_most and s.surface_heat == "temp" and deck.get_str("erf.most.pot_surf_temp") != rows[0][1]:
            failures.append("%s: erf.most.pot_surf_temp is not written as the sounding's surface theta" % label)
    if s.zlo_most and s.roughness == "const" and s.surface_heat != "file":
        close("erf.most.z0", s.z0)
    if s.zlo_most and s.surface_heat == "flux":
        close("erf.most.surf_temp_flux", s.heat_flux)
    if s.zlo_most and s.surface_heat == "file" and deck.get_int("erf.is_land") != 0:
        failures.append("%s: the sea-surface table needs erf.is_land = 0" % label)
    if s.coriolis:
        close("erf.latitude", s.latitude)
    if s.setup == "bubble" and s.dycore == "anelastic" and deck.get_str("prob.perturb_rho") != "false":
        failures.append("%s: an anelastic bubble needs prob.perturb_rho = false" % label)
    for key in ("erf.rayleigh_dampcoef", "erf.sponge_strength"):
        rate = deck.get_real(key)
        if rate is not None and rate * s.dt_est > 0.25:
            failures.append("%s: %s x time step is %g, above 0.25" % (label, key, rate * s.dt_est))
    if deck.get_ints("amr.n_cell") != list(s.n):
        failures.append("%s: amr.n_cell does not hold the requested cells" % label)
    closures = sum(1 for k in ("erf.les_type", "erf.rans_type", "erf.pbl_type")
                   if (deck.get_str(k) or "None") != "None")
    expected = 0 if s.turbulence in ("none", "dns") else 1
    if closures != expected:
        failures.append("%s: %d turbulence closures are active, expected %d" % (label, closures, expected))
    if ("plt" in s.outputs) != deck.has("erf.plot_vars_1"):
        failures.append("%s: erf.plot_vars_1 does not follow --output" % label)


# cadence keys set to every step in the ERF runs, so each active output fires after step 0
_EVERY_STEP = ("erf.plot_int_1", "erf.plot2d_int_1", "erf.profile_int", "erf.sum_interval", "erf.line_sampling_interval",
               "erf.plane_sampling_interval", "erf.station_sampling_interval", "erf.subvol_int")


def _run_erf(exe: str, s, deck, outdir: str, failures, label) -> None:
    cmd = [os.path.abspath(exe), s.name, "max_step=3", "amrex.fpe_trap_invalid=1", "erf.check_int=-1"]
    cmd += ["%s=1" % key for key in _EVERY_STEP if deck.is_set(key)]
    try:
        res = subprocess.run(cmd, cwd=outdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True, timeout=600)
    except (OSError, subprocess.TimeoutExpired) as exc:
        failures.append("%s: ERF did not run: %s" % (label, exc))
        return
    if res.returncode != 0 or "finalized" not in res.stdout:
        lines = res.stdout.strip().splitlines()
        tail = " | ".join([ln for ln in lines if "Abort" in ln or "Assert" in ln or "Error" in ln][:2] or lines[-4:])
        failures.append("%s: ERF stopped (exit %d): %s" % (label, res.returncode, tail))
        return
    # AMReX lists, at the end of a run, every input that nothing read
    if "Unused ParmParse Variables" in res.stdout:
        unused = res.stdout.split("Unused ParmParse Variables:", 1)[1].split("AMReX (", 1)[0].split()
        failures.append("%s: ERF did not read %s" % (label, " ".join(u for u in unused if "::" in u and "[" in u)))


def main(schema_path: Optional[str] = None) -> int:
    failures: List[str] = []
    schema = None
    if schema_path:
        with open(schema_path) as fh:
            schema = json.load(fh)
    read_names = _source_literals()
    if not read_names:
        failures.append("the key names read under Source/ were not found; run from an ERF source tree")
    from erf_input_checker.checks.keys import documented_keys
    documented = documented_keys()
    if not documented:
        failures.append("Docs/sphinx_doc/Inputs.rst was not found; run from an ERF source tree")
    exe = os.environ.get("ERF_CREATEINPUT_EXE")
    ran = 0
    with tempfile.TemporaryDirectory() as top:
        for n, options in enumerate(CASES):
            label = " ".join(options) or "(defaults)"
            outdir = os.path.join(top, "case%03d" % n)
            # small grids keep the optional ERF runs short; a case that sets its own grid keeps it
            # (the density current and the hill keep cubic cells: the multigrid Poisson solve of the
            # anelastic cases did not converge on 1600 x 100 x 400 m and 250 x 25 x 125 m cells)
            setup = parse(options).setup
            if "--nx" in options or "--domain" in options:
                small = []
            elif setup == "hill":
                small = ["--domain", "1600", "400", "800", "--nx", "32", "--ny", "8", "--nz", "16"]
            elif setup == "density_current":
                small = ["--domain", "3200", "100", "1600", "--nx", "32", "--ny", "1", "--nz", "16"]
            else:
                small = ["--nx", "16", "--ny", "16", "--nz", "16"]
            try:
                s = _spec(options + small, outdir)
                deck, files = build(s)
                write(s, deck, files, case_options(parse(options + small)))
            except (CaseError, ValueError) as exc:
                failures.append("%s: not written: %s" % (label, exc))
                continue
            _check_keys(deck, schema, read_names, documented, failures, label)
            # the case file written beside the deck must give the same deck again
            again = resolve(parse(["--case-file", os.path.join(outdir, "case.json"), "--outdir", outdir,
                                   "--name", s.name]))
            if build(again)[0].render() != deck.render():
                failures.append("%s: case.json does not reproduce the deck" % label)
            _check_aux(s, outdir, failures, label)
            _check_deck(s, outdir, failures, label, needed=[f.key for f in files if not f.generated])
            if sorted(os.listdir(outdir)) != sorted([s.name] + [f.name for f in files if f.generated]
                                                    + (["FILES_TO_REPLACE.txt"] if files else []) + ["case.json"]):
                failures.append("%s: the files written do not match the files the deck names" % label)
            # a case that waits for files the tool cannot generate cannot start
            # (NetCDF plotfiles need a NetCDF build: run that case only when ERF_CREATEINPUT_NETCDF is set)
            # ERF_CREATEINPUT_FILTER runs only the cases whose options hold that text (a build with
            # wind farms, for example, runs only wind-farm cases)
            if exe and ("netcdf" not in s.outputs or os.environ.get("ERF_CREATEINPUT_NETCDF")) \
                    and all(f.generated for f in files) and os.environ.get("ERF_CREATEINPUT_FILTER", "") in label:
                _run_erf(exe, s, deck, outdir, failures, label)
                ran += 1
        # a second write into the same directory must be refused without --force
        s = _spec(["--setup", "abl"], os.path.join(top, "case000"))
        try:
            write(s, *build(s))
            failures.append("an existing case was overwritten without --force")
        except CaseError:
            pass
        # an option on the command line replaces the preset's value, and the preset fills the rest
        s = _spec(["--preset", "abl_convective_les", "--turbulence", "smagorinsky"], os.path.join(top, "over"))
        if s.turbulence != "smagorinsky" or s.surface_heat != "flux" or not s.coriolis:
            failures.append("--preset with --turbulence: the command line must win and the preset fill the rest")
        with open(os.path.join(top, "c.json"), "w") as fh:
            json.dump({"setup": "hill", "dycore": "anelastic", "z0": 0.05}, fh)
        s = _spec(["--preset", "abl_rans", "--case-file", os.path.join(top, "c.json"), "--z0", "0.2"],
                  os.path.join(top, "over"))
        if (s.setup, s.dycore, s.z0, s.turbulence) != ("hill", "anelastic", 0.2, "keqn"):
            failures.append("--case-file must replace preset values, and the command line both")
        for bad in ({"nosuch": 1}, {"setup": "moon"}):
            with open(os.path.join(top, "bad.json"), "w") as fh:
                json.dump(bad, fh)
            try:
                parse(["--case-file", os.path.join(top, "bad.json")])
                failures.append("a case file holding %s was accepted" % bad)
            except OptionError:
                pass
        # an active key written twice must be refused: ParmParse would take the last one silently
        twice = DeckText()
        sec = twice.section("t")
        sec.set("erf.cfl", 0.5)
        try:
            sec.set("erf.cfl", 0.8)
            failures.append("a key was set twice in one deck without an error")
        except ValueError:
            pass
        for options, needle in REFUSED:
            label = " ".join(options)
            try:
                build(_spec(options, os.path.join(top, "refused")))
                failures.append("%s: was not refused" % label)
            except CaseError as exc:
                if needle not in str(exc):
                    failures.append("%s: the message does not name %s: %s" % (label, needle, exc))
    for line in failures:
        print("FAIL " + line)
    print("erfcreateinput self-test: %d cases written, %d refused, %d run in ERF, %d failure(s)"
          % (len(CASES), len(REFUSED), ran, len(failures)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
