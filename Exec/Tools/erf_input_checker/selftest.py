"""Built-in tests: each check is proven against a deck that trips it and a
sibling deck that does not (``--self-test``).

Every case writes a small deck to a temporary directory, runs the checkers
and asserts on the (severity, key) pairs in the report.  A check whose
"bad" deck does not produce the expected finding, or whose "good" deck does,
fails the test, so a check cannot be silently broken by a later edit.
"""

import os
import shutil
import tempfile
from typing import List, Tuple

from .checker import run_all_checks
from .context import Options, RunContext
from .deck import parse_deck
from .findings import Report, Severity

BASE = """
amr.n_cell = 64 64 32
amr.max_level = 0
geometry.prob_lo = 0 0 0
geometry.prob_hi = 6400 6400 2000
geometry.is_periodic = 1 1 0
zlo.type = SlipWall
zhi.type = SlipWall
erf.init_type = Uniform
erf.fixed_dt = 1.0
max_step = 10
erf.plot_int_1 = 10
erf.plot_vars_1 = density theta
"""

# (name, extra deck text, expected (severity, key) present, expected absent)
CASES: List[Tuple[str, str, List[Tuple[str, str]], List[Tuple[str, str]]]] = [
    ("clean", "", [], [("ERROR", "*"), ("WARNING", "*")]),
    ("odd max_grid_size", "amr.max_grid_size = 15", [("ADVICE", "amr.max_grid_size")], [("ERROR", "amr.max_grid_size")]),
    ("odd mgs in z is fine", "amr.max_grid_size_z = 15", [], [("ADVICE", "amr.max_grid_size")]),
    ("even mgs", "amr.max_grid_size = 16", [], [("ERROR", "amr.max_grid_size")]),
    ("uneven boxes", "amr.max_grid_size = 48", [("ADVICE", "amr.max_grid_size")], []),
    ("blocking factor not pow2", "amr.blocking_factor = 6", [("ERROR", "amr.blocking_factor")], []),
    ("n_cell not multiple of bf", "amr.blocking_factor = 16\namr.n_cell = 72 64 32", [("ERROR", "amr.n_cell")], []),
    ("bf in z skipped", "amr.blocking_factor = 16\namr.n_cell = 64 64 36", [], [("ERROR", "amr.n_cell")]),
    ("odd n_cell", "amr.n_cell = 63 64 32", [("ADVICE", "amr.n_cell")], []),
    ("ranks imbalance", "amr.max_grid_size = 32", [("ADVICE", "amr.max_grid_size")], [], ),
    ("cfl ignored with fixed_dt", "erf.cfl = 0.5", [("WARNING", "erf.cfl")], []),
    ("cfl above one is a warning", "erf.cfl = 1.5\nerf.fixed_dt = -1", [("WARNING", "erf.cfl")], [("ERROR", "erf.cfl")]),
    ("cfl zero with fixed_dt is fine", "erf.cfl = 0", [], [("ERROR", "erf.cfl")]),
    ("cfl zero without fixed_dt", "erf.cfl = 0\nerf.fixed_dt = -1", [("ERROR", "erf.cfl")], []),
    ("cfl zero no substep no crash", "erf.cfl = 0\nerf.fixed_dt = -1\nerf.substepping_type = None", [("ERROR", "erf.cfl")], [("ERROR", "internal")]),
    ("substepping_cfl zero", "erf.substepping_cfl = 0", [("ERROR", "erf.substepping_cfl")], [("ERROR", "internal")]),
    ("mri ratio zero is unset", "erf.fixed_mri_dt_ratio = 0", [], [("ERROR", "erf.fixed_mri_dt_ratio")]),
    ("max_level not an integer", "amr.max_level = two", [("ERROR", "amr.max_level")], [("ERROR", "internal")]),
    ("prob_hi and prob_extent", "geometry.prob_extent = 6400 6400 2000", [("ERROR", "geometry.prob_extent")], []),
    ("check_int and check_per", "amr.check_int = 10\namr.check_per = 1.0", [("WARNING", "amr.check_int")], [("ERROR", "amr.check_int")]),
    ("key case", "erf.Fixed_dt = 2.0", [("WARNING", "erf.Fixed_dt")], []),
    ("wrong prefix", "amr.fixed_dt = 2.0", [("WARNING", "amr.fixed_dt")], [("INFO", "amr.fixed_dt")]),
    ("molec_diff_type case", "erf.molec_diff_type = constant", [("ERROR", "erf.molec_diff_type")], []),
    ("advection case", "erf.dycore_horiz_adv_type = upwind_5th", [("ERROR", "erf.dycore_horiz_adv_type")], []),
    ("les case", "erf.les_type = smagorinsky", [("ERROR", "erf.les_type")], []),
    ("pbl case insensitive", "erf.pbl_type = mynn25\nzlo.type = surface_layer", [], [("ERROR", "erf.pbl_type")]),
    ("moisture WSM6", "erf.moisture_model = WSM6\nerf.use_gravity = true", [], [("ERROR", "erf.moisture_model"), ("ADVICE", "erf.moisture_model")]),
    ("moisture unknown", "erf.moisture_model = Thompson", [("ERROR", "erf.moisture_model")], []),
    ("use_real_bcs on ideal is a warning", "erf.use_real_bcs = true", [("WARNING", "erf.use_real_bcs")], [("ERROR", "erf.use_real_bcs")]),
    ("wrfinput flat terrain", "erf.init_type = WRFInput\nerf.nc_init_file_0 = wrfinput_d01\nerf.use_gravity = true\nerf.moisture_model = Kessler\ngeometry.is_periodic = 0 0 0",
     [("ERROR", "erf.terrain_type")], []),
    ("wrfinput dry", "erf.init_type = WRFInput\nerf.nc_init_file_0 = wrfinput_d01\nerf.use_gravity = true\nerf.terrain_type = StaticFittedMesh\ngeometry.is_periodic = 0 0 0",
     [("ERROR", "erf.moisture_model")], []),
    ("odd mri ratio", "erf.fixed_mri_dt_ratio = 3", [("ERROR", "erf.fixed_mri_dt_ratio")], []),
    ("even mri ratio", "erf.fixed_mri_dt_ratio = 4", [], [("ERROR", "erf.fixed_mri_dt_ratio")]),
    ("fast dt non-integer ratio", "erf.fixed_fast_dt = 0.3", [("ERROR", "erf.fixed_fast_dt")], []),
    ("fast dt odd integer ratio", "erf.fixed_fast_dt = 0.2", [("ADVICE", "erf.fixed_fast_dt")], [("ERROR", "erf.fixed_fast_dt")]),
    ("fast dt even ratio", "erf.fixed_fast_dt = 0.125", [], [("ERROR", "erf.fixed_fast_dt"), ("ADVICE", "erf.fixed_fast_dt")]),
    ("over-specified dt", "erf.fixed_fast_dt = 0.25\nerf.fixed_mri_dt_ratio = 6", [("ERROR", "erf.fixed_mri_dt_ratio")], []),
    ("too many substeps", "erf.fixed_dt = 20", [("WARNING", "erf.fixed_dt")], []),
    ("advective cfl too big", "erf.fixed_dt = 6", [("WARNING", "erf.fixed_dt")], []),
    ("no substep acoustic cfl", "erf.substepping_type = None\nerf.fixed_dt = 1.0", [("WARNING", "erf.fixed_dt")], []),
    ("no substep ok", "erf.substepping_type = None\nerf.fixed_dt = 0.1", [], [("WARNING", "erf.fixed_dt", "timestep.compressible")]),
    ("fast keys with no substepping", "erf.substepping_type = None\nerf.fixed_mri_dt_ratio = 4", [("WARNING", "erf.fixed_mri_dt_ratio")], []),
    ("anelastic ignores fast keys", "erf.anelastic = 1\nerf.use_fft = true\nerf.fixed_mri_dt_ratio = 4", [("WARNING", "erf.fixed_mri_dt_ratio")], []),
    ("anelastic without fft", "erf.anelastic = 1", [("ADVICE", "erf.use_fft")], []),
    ("anelastic with fft", "erf.anelastic = 1\nerf.use_fft = true", [], [("ADVICE", "erf.use_fft")]),
    ("anelastic advective cfl", "erf.anelastic = 1\nerf.use_fft = true\nerf.fixed_dt = 6", [("WARNING", "erf.fixed_dt")], []),
    ("diffusive limit with fixed dt", "erf.diffusive_dt_limit = true\nerf.les_type = Smagorinsky", [("ERROR", "erf.diffusive_dt_limit")], []),
    ("diffusive cfl alone", "erf.diffusive_cfl = 0.5", [("ERROR", "erf.diffusive_cfl")], []),
    ("never ends", "max_step = -1", [("ERROR", "max_step")], []),
    ("stop_time ok", "max_step = -1\nstop_time = 100", [], [("ERROR", "max_step")]),
    ("missing face type", "geometry.is_periodic = 0 1 0\nxlo.type = inflow\nxlo.velocity = 1 0 0", [("ERROR", "xhi.type")], []),
    ("face keys on periodic dir", "xlo.type = inflow", [("WARNING", "xlo.type")], []),
    ("inflow without velocity", "geometry.is_periodic = 0 1 0\nxlo.type = inflow\nxhi.type = outflow", [("ERROR", "xlo.type")], []),
    ("bad face type", "zhi.type = lid", [("ERROR", "zhi.type")], []),
    ("velocity on outflow", "geometry.is_periodic = 0 1 0\nxlo.type = inflow\nxlo.velocity = 1 0 0\nxhi.type = outflow\nxhi.velocity = 1 0 0", [("WARNING", "xhi.velocity")], []),
    ("most keys without surface layer", "erf.most.z0 = 0.1", [("WARNING", "erf.most.z0")], []),
    ("pbl without surface layer", "erf.pbl_type = MYNN25", [("ERROR", "erf.pbl_type")], []),
    ("surface layer no diffusion", "zlo.type = surface_layer", [("WARNING", "zlo.type")], []),
    ("z0 above zref", "zlo.type = surface_layer\nerf.les_type = Smagorinsky\nerf.most.z0 = 40", [("ERROR", "erf.most.z0")], []),
    ("retired key", "erf.use_terrain = true", [("ERROR", "erf.use_terrain")], []),
    ("retired value", "erf.init_type = Real", [("ERROR", "erf.init_type")], []),
    ("misspelt key", "erf.dynamicViscosity = 1e-5", [("WARNING", "erf.dynamicViscosity")], []),
    ("duplicate differing", "erf.fixed_dt = 2.0", [("WARNING", "erf.fixed_dt")], []),
    ("gravity missing", "erf.init_type = Input_Sounding", [("WARNING", "erf.use_gravity")], []),
    ("real keys in ideal deck", "erf.nc_bdy_file = wrfbdy_d01", [("WARNING", "erf.nc_bdy_file")], []),
    ("real init missing file", "erf.init_type = WRFInput\nerf.use_gravity = true\nerf.terrain_type = StaticFittedMesh\nerf.moisture_model = Kessler",
     [("ERROR", "erf.nc_init_file_0")], []),
    ("real init file absent", "erf.init_type = WRFInput\nerf.nc_init_file_0 = wrfinput_d01\nerf.nc_bdy_file = wrfbdy_d01\nerf.use_gravity = true\nerf.terrain_type = StaticFittedMesh\nerf.moisture_model = Kessler\ngeometry.is_periodic = 0 0 0",
     [("ERROR", "erf.nc_init_file_0"), ("ERROR", "erf.nc_bdy_file")], [("ERROR", "xlo.type")]),
    ("periodic with real bcs", "erf.init_type = WRFInput\nerf.nc_init_file_0 = wrfinput_d01\nerf.nc_bdy_file = wrfbdy_d01\nerf.use_gravity = true\nerf.terrain_type = StaticFittedMesh\nerf.moisture_model = Kessler",
     [("ERROR", "geometry.is_periodic")], []),
    ("les and rans together", "erf.les_type = Smagorinsky\nerf.rans_type = kEqn", [("WARNING", "erf.les_type")], []),
    ("Cs without smagorinsky", "erf.Cs = 0.1", [("WARNING", "erf.Cs")], []),
    ("centered without num diff", "erf.dycore_horiz_adv_type = Centered_2nd", [("ADVICE", "erf.dycore_horiz_adv_type")], []),
    ("centered with num diff", "erf.dycore_horiz_adv_type = Centered_2nd\nerf.num_diff_coeff = 0.05", [], [("ADVICE", "erf.dycore_horiz_adv_type")]),
    ("moist adv without moisture", "erf.moistscal_horiz_adv_type = WENO5", [("WARNING", "erf.moistscal_horiz_adv_type")], []),
    ("geostrophic without coriolis", "erf.abl_driver_type = GeostrophicWind\nerf.abl_geo_wind = 10 0 0", [("WARNING", "erf.abl_driver_type")], []),
    ("coriolis default latitude", "erf.use_coriolis = true", [("WARNING", "erf.latitude")], []),
    ("rayleigh depth", "erf.rayleigh_damp_W = true\nerf.rayleigh_zdamp = 3000", [("WARNING", "erf.rayleigh_zdamp")], []),
    ("stretching without initial_dz", "erf.grid_stretching_ratio = 1.05", [("ERROR", "erf.initial_dz")], []),
    ("stretched top mismatch", "erf.grid_stretching_ratio = 1.05\nerf.initial_dz = 10", [("WARNING", "erf.initial_dz")], []),
    ("stretched top matches", "erf.grid_stretching_ratio = 1.05\nerf.initial_dz = %.10g" % (2000 * 0.05 / (1.05 ** 32 - 1)), [], [("WARNING", "erf.initial_dz"), ("ADVICE", "erf.initial_dz")]),
    ("z levels count", "erf.terrain_z_levels = 0 100 200", [("ERROR", "erf.terrain_z_levels")], []),
    ("SAM with stretching", "erf.moisture_model = SAM\nerf.grid_stretching_ratio = 1.05\nerf.initial_dz = 10", [("ERROR", "erf.moisture_model")], []),
    ("refinement without indicator", "amr.max_level = 1", [("WARNING", "erf.refinement_indicators")], []),
    ("indicator without keys", "amr.max_level = 1\nerf.refinement_indicators = box1", [("ERROR", "erf.refinement_indicators")], []),
    ("box outside domain", "amr.max_level = 1\nerf.refinement_indicators = box1\nerf.box1.in_box_lo = 0 0 0\nerf.box1.in_box_hi = 9000 100 100", [("WARNING", "erf.box1.in_box_lo")], []),
    ("ref ratio 1", "amr.max_level = 1\namr.ref_ratio = 1\nerf.refinement_indicators = box1\nerf.box1.in_box_lo = 0 0 0\nerf.box1.in_box_hi = 100 100 100", [("ERROR", "amr.ref_ratio")], []),
    ("two-way with cf_width", "amr.max_level = 1\nerf.cf_width = 4\nerf.refinement_indicators = box1\nerf.box1.in_box_lo = 0 0 0\nerf.box1.in_box_hi = 100 100 100", [("ERROR", "erf.cf_width")], []),
    ("plot int and per", "erf.plot_per_1 = 1.0", [("ERROR", "erf.plot_int_1")], []),
    ("plot vars without int", "erf.plot_vars_2 = density", [("WARNING", "erf.plot_vars_2")], []),
    ("legacy plot key", "erf.plot_int = 5", [("WARNING", "erf.plot_int")], []),
    ("no checkpoints long run", "max_step = 50000", [("ADVICE", "amr.check_int")], []),
    ("restart missing", "amr.restart = chk00100", [("WARNING", "amr.restart")], []),
    ("vertical split with pbl", "amr.no_box_split_dir = -1\nerf.substepping_type = None\namr.max_grid_size = 16\nzlo.type = surface_layer\nerf.pbl_type = MYNN25",
     [("ERROR", "amr.max_grid_size")], []),
    ("implicit substep needs nbsd 2", "amr.no_box_split_dir = -1", [("ERROR", "amr.no_box_split_dir")], []),
    ("terrain file missing", "erf.terrain_type = StaticFittedMesh\nerf.terrain_file_name = hill.txt", [("ERROR", "erf.terrain_file_name")], []),
    ("terrain source with flat mesh", "erf.terrain_file_name = hill.txt", [("WARNING", "erf.terrain_file_name")], []),
    ("aspect ratio", "erf.les_type = Smagorinsky\namr.n_cell = 64 64 200", [("ADVICE", "erf.les_type")], []),
]

# file-backed cases: (name, deck text, files {name: content}, expected present, expected absent)
FILE_CASES = [
    ("terrain file short in x",
     BASE + "erf.terrain_type = StaticFittedMesh\nerf.terrain_file_name = hill.txt\n",
     {"hill.txt": "3 2\n0 1000 2000\n0 6400\n" + " ".join(["10"] * 6) + "\n"},
     [("ERROR", "erf.terrain_file_name")], []),
    ("terrain file covers domain",
     BASE + "erf.terrain_type = StaticFittedMesh\nerf.terrain_file_name = hill.txt\n",
     {"hill.txt": "3 2\n0 3200 6400\n0 6400\n" + " ".join(["10"] * 6) + "\n"},
     [("INFO", "erf.terrain_file_name")], [("ERROR", "erf.terrain_file_name")]),
    ("terrain file wrong count",
     BASE + "erf.terrain_type = StaticFittedMesh\nerf.terrain_file_name = hill.txt\n",
     {"hill.txt": "3 2\n0 3200 6400\n0 6400\n10 10\n"},
     [("ERROR", "erf.terrain_file_name")], []),
    ("sounding surface theta mismatch",
     BASE + "erf.init_type = Input_Sounding\nerf.use_gravity = true\nerf.input_sounding_file = snd\nzlo.type = surface_layer\nerf.les_type = Smagorinsky\nerf.most.surf_temp = 305\n",
     {"snd": "1000 300 0\n0 300 0 5 0\n1000 302 0 5 0\n3000 310 0 5 0\n"},
     [("ERROR", "erf.most.surf_temp")], []),
    ("sounding below domain top",
     BASE + "erf.init_type = Input_Sounding\nerf.use_gravity = true\nerf.input_sounding_file = snd\n",
     {"snd": "1000 300 0\n0 300 0 5 0\n1000 302 0 5 0\n"},
     [("WARNING", "erf.input_sounding_file")], []),
    ("sounding fine",
     BASE + "erf.init_type = Input_Sounding\nerf.use_gravity = true\nerf.input_sounding_file = snd\n",
     {"snd": "1000 300 0\n0 300 0 5 0\n3000 310 0 5 0\n"},
     [], [("WARNING", "erf.input_sounding_file"), ("ERROR", "erf.input_sounding_file")]),
    ("sounding with two columns on the first line",
     BASE + "erf.init_type = Input_Sounding\nerf.use_gravity = true\nerf.input_sounding_file = snd\n",
     {"snd": "1000 300\n0 300 0 5 0\n3000 310 0 5 0\n"},
     [("ERROR", "erf.input_sounding_file")], [("ERROR", "internal")]),
    ("every file of a multi-file key is checked",
     BASE + "erf.init_type = Input_Sounding\nerf.use_gravity = true\nerf.input_sounding_file = snd missing_snd\nerf.input_sounding_time = 0 3600\nerf.nudging_from_input_sounding = true\n",
     {"snd": "1000 300 0\n0 300 0 5 0\n3000 310 0 5 0\n"},
     [("ERROR", "erf.input_sounding_file")], []),
    ("spacedim branch not taken",
     BASE + "#if AMREX_SPACEDIM == 2\nerf.fixed_dt = 99\n#else\nerf.cfl = 0.5\nerf.fixed_dt = -1\n#endif\n",
     {},
     [("INFO", "erf.cfl")], [("WARNING", "erf.fixed_dt", "timestep.compressible"), ("WARNING", "erf.cfl")]),
    ("namelist skipped",
     BASE + "&fortran_block\n  nx = 10,\n/\n",
     {},
     [], [("ERROR", "*"), ("WARNING", "*")]),
    ("two definitions on one line",
     BASE + "erf.plot_int_2 = 5 erf.plot_vars_2 = density\n",
     {},
     [], [("WARNING", "erf.plot_vars_2"), ("ERROR", "*")]),
    ("include found beside the deck only",
     "FILE = sub/base\nerf.fixed_dt = 0.5\n",
     {"sub/base": BASE},
     [("WARNING", "FILE")], []),
    ("include override is info",
     "FILE = base\nerf.fixed_dt = 0.5\n",
     {"base": BASE},
     [("INFO", "erf.fixed_dt")], [("WARNING", "erf.fixed_dt")]),
]


def _run(path: str) -> Report:
    deck = parse_deck(path)
    rep = Report(path)
    run_all_checks(deck, RunContext(deck, Options(ranks=3)), rep)
    return rep


def _has(rep: Report, sev: str, key: str, checker: str = None) -> bool:
    for f in rep.findings:
        if f.severity.name == sev and (key == "*" or f.key == key) and (checker is None or f.checker == checker):
            return True
    return False


def main() -> int:
    tmp = tempfile.mkdtemp(prefix="erf_input_checker_")
    failures = []
    try:
        for name, extra, present, absent in CASES:
            path = os.path.join(tmp, "inputs")
            with open(path, "w") as fh:
                fh.write(BASE + extra + "\n")
            _check_case(name, path, present, absent, failures)
        for name, text, files, present, absent in FILE_CASES:
            case_dir = os.path.join(tmp, "case")
            shutil.rmtree(case_dir, ignore_errors=True)
            os.makedirs(case_dir)
            for fn, content in files.items():
                full = os.path.join(case_dir, fn)
                os.makedirs(os.path.dirname(full), exist_ok=True)
                with open(full, "w") as fh:
                    fh.write(content)
            path = os.path.join(case_dir, "inputs")
            with open(path, "w") as fh:
                fh.write(text)
            _check_case(name, path, present, absent, failures)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    total = len(CASES) + len(FILE_CASES)
    if failures:
        print("\n".join(failures))
        print("self-test: %d of %d cases FAILED" % (len(failures), total))
        return 1
    print("self-test: %d cases passed" % total)
    return 0


def _check_case(name, path, present, absent, failures) -> None:
    rep = _run(path)
    for item in present:
        if not _has(rep, *item):
            failures.append("FAIL %-32s expected %s; got:\n%s" % (name, item, rep.format(Severity.INFO)))
    for item in absent:
        if _has(rep, *item):
            failures.append("FAIL %-32s did not expect %s; got:\n%s" % (name, item, rep.format(Severity.INFO)))
