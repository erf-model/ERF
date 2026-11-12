"""Option ranges and the combinations ERF refuses at start-up.

Each rule returns a message that names the options involved, so the case is
refused before a deck is written.  Nothing is clamped or changed silently.
The ERF source line behind each combination rule is cited beside it.
"""

import math
import os
from typing import Callable, List, Tuple

from .options import OUTPUTS, PBL, REAL_SETUPS
from .recipes.models import ICE_MODELS, RAIN_MODELS

Rule = Tuple[Callable, str]

# names the tool writes beside the deck
RESERVED_NAMES = ("FILES_TO_REPLACE.txt", "case.json", "input_sounding", "roughness.txt", "sfc.txt",
                  "geo_wind_table.txt", "terrain.txt", "terrain_usgs.txt", "buildings.txt", "forest.txt",
                  "windturbines_loc_x_y.txt", "windturbines_spec.tbl")


def _ranges(s) -> List[str]:
    bad = []
    numbers = [("--dt", s.dt), ("--stop-time", s.stop_time), ("--theta", s.theta), ("--z0", s.z0),
               ("--latitude", s.latitude), ("--heat-flux", s.heat_flux), ("--dz0", s.dz0),
               ("--hill-height", s.hill_height), ("--hill-width", s.hill_width),
               ("--wind", s.wind[0]), ("--wind", s.wind[1])]
    numbers += [("--domain", v) for v in (s.domain or [])]
    for name, value in numbers:
        if value is not None and not math.isfinite(value):
            bad.append("%s = %s: give a finite number" % (name, value))
    if bad:
        return bad
    for name, value, least in (("--nx", s.n[0], 1), ("--ny", s.n[1], 1), ("--nz", s.n[2], 2)):
        if value < least:
            bad.append("%s = %d: the cell count must be at least %d" % (name, value, least))
    if s.domain is not None:
        for axis, ext in zip("xyz", s.domain):
            if not ext > 0.0:
                bad.append("--domain: the %s extent %g m must be positive" % (axis, ext))
    if s.dt is not None and not s.dt > 0.0:
        bad.append("--dt = %g: the time step must be positive" % s.dt)
    if not s.stop_time > 0.0:
        bad.append("--stop-time = %g: the run length must be positive" % s.stop_time)
    if not s.theta > 0.0:
        bad.append("--theta = %g: the potential temperature must be positive" % s.theta)
    if not s.z0 > 0.0:
        bad.append("--z0 = %g: the roughness length must be positive" % s.z0)
    if not -90.0 <= s.latitude <= 90.0:
        bad.append("--latitude = %g: must be between -90 and 90 degrees" % s.latitude)
    if not s.hill_height > 0.0:
        bad.append("--hill-height = %g: the hill height must be positive" % s.hill_height)
    if not s.hill_width > 0.0:
        bad.append("--hill-width = %g: the hill half-width must be positive" % s.hill_width)
    if s.terrain == "analytic" and s.hill_height > 0.5 * (s.hi[2] - s.lo[2]):
        bad.append("--hill-height = %g: must be at most half the domain depth, %g m"
                   % (s.hill_height, 0.5 * (s.hi[2] - s.lo[2])))
    if s.dz0 is not None and not 0.0 < s.dz0 < s.dx[2]:
        bad.append("--dz0 = %g: the first-cell height must be positive and below the uniform cell height %g m"
                   % (s.dz0, s.dx[2]))
    elif s.dz0 is not None and s.stretch is None:
        # Source/Utils/ERF_InitZLevels.cpp:73 asserts when the Newton iteration of geometric_stretch fails
        bad.append("--dz0 = %g: too thin for %d cells over %g m; ERF's solve for the stretching ratio does not "
                   "converge. Use a thicker first cell or more cells" % (s.dz0, s.n[2], s.hi[2] - s.lo[2]))
    reserved = RESERVED_NAMES
    if s.name in reserved or s.name in (".", "..") or os.sep in s.name:
        bad.append("--name = '%s': the name is used for a generated file or is not a file name" % s.name)
    if not 0 <= s.amr <= 2:
        bad.append("--amr = %d: the number of refined levels must be 0, 1 or 2" % s.amr)
    bad += s.refine_errors
    if s.refine_box is not None:
        box = s.refine_box
        if len(box) not in (4, 6) or not all(math.isfinite(v) for v in box):
            bad.append("--refine-box: give XLO YLO XHI YHI or XLO YLO ZLO XHI YHI ZHI, as finite numbers")
        elif any(box[k] >= box[k + len(box) // 2] for k in range(len(box) // 2)):
            bad.append("--refine-box: each upper value must exceed the lower one")
    for name in s.outputs:
        if name not in OUTPUTS:
            bad.append("--output: '%s' is not one of %s" % (name, ", ".join(OUTPUTS)))
    if not s.name or "/" in s.name:
        bad.append("--name = '%s': give a file name without a directory" % s.name)
    return bad


RULES: List[Rule] = [
    # Source/ERF.cpp:1903-1911: a surface_layer face needs a diffusive closure
    (lambda s: s.surface == "most" and s.turbulence == "none",
     "--surface most needs a turbulence closure; --turbulence none has no diffusion and ERF stops"),
    # Source/DataStructs/ERF_TurbStruct.H:385-391: these schemes need zlo.type = surface_layer
    (lambda s: s.turbulence in PBL and s.turbulence != "myj" and s.surface != "most",
     "--turbulence {turbulence} needs --surface most (the scheme reads the surface-layer fluxes)"),
    # Source/PBL/ERF_ComputeDiffusivityYSU.cpp:131-133: YSU stops unless the surface layer is stable
    (lambda s: s.turbulence == "ysu" and not (s.surface_heat == "flux" and s.heat_flux < 0.0),
     "--turbulence ysu supports stable conditions only: add --surface-heat flux with a negative "
     "--heat-flux, or use --turbulence ysunew"),
    # ERF_DataStruct.H:1739-1742: abl_geo_wind is read only when use_coriolis is true
    (lambda s: s.forcing in ("geostrophic", "geo_table") and not s.coriolis,
     "--forcing {forcing} needs --coriolis; without it ERF never reads the geostrophic wind"),
    (lambda s: s.roughness == "file" and s.surface != "most",
     "--roughness file needs --surface most"),
    (lambda s: s.surface_heat != "adiabatic" and s.surface != "most",
     "--surface-heat {surface_heat} needs --surface most"),
    # Source/Initialization/ERF_Init1D.cpp:74: a uniform base state has no hydrostatic balance,
    # so gravity is off and a surface heat source would have no buoyancy to act through
    (lambda s: s.surface_heat != "adiabatic" and s.init == "uniform",
     "--surface-heat {surface_heat} needs --init sounding or isentropic; --init uniform runs "
     "without gravity, so the heating would not drive any buoyancy"),
    # both setups are driven by buoyancy, which a uniform base state without gravity does not have
    (lambda s: s.setup in ("bubble", "density_current") and s.init == "uniform",
     "--setup {setup} needs --init isentropic or sounding; --init uniform runs without gravity"),
    # the roughness file has one row per node of the level-0 grid
    (lambda s: s.roughness == "file" and (s.n[0] < 2 or s.n[1] < 2),
     "--roughness file needs at least 2 cells in x and y"),
    (lambda s: s.terrain == "none" and s.terrain_method != "fitted",
     "--terrain-method {terrain_method} needs terrain; add --terrain analytic, file or file_usgs"),
    # ERF_DataStruct.H ~2180: Deardorff and kEqn stop with terrain_type = EB
    (lambda s: s.method == "eb" and s.turbulence in ("deardorff", "keqn"),
     "--terrain-method eb does not support --turbulence {turbulence}"),
    # ERF_SurfaceLayer.H:612-614, 747-755: no roughness file or surface table with EB
    (lambda s: s.method == "eb" and s.surface == "most" and (s.roughness == "file" or s.surface_heat == "file"),
     "--terrain-method eb does not support --roughness file or --surface-heat file"),
    # ERF_DataStruct.H:1053-1055
    (lambda s: s.method == "eb" and s.buildings_on,
     "--buildings file uses immersed forcing, which ERF does not combine with --terrain-method eb"),
    # stretching needs the terrain-fitted (variable dz) mesh; immersed forcing and EB keep constant dz
    (lambda s: s.stretched and s.method in ("immersed", "eb"),
     "--dz0 needs --terrain-method fitted (or no terrain); {terrain_method} keeps uniform cells"),
    (lambda s: s.stretched and s.buildings_on,
     "--dz0 cannot be combined with --buildings file; immersed forcing keeps uniform cells"),
    # the immersed surface takes one roughness length (erf.if_z0) and a flux or a temperature
    (lambda s: s.method == "immersed" and s.surface == "most" and (s.roughness == "file" or s.surface_heat == "file"),
     "--terrain-method immersed takes one roughness length and a flux or temperature; "
     "--roughness file and --surface-heat file are not available"),
    (lambda s: s.method == "immersed" and s.turbulence in PBL,
     "--turbulence {turbulence} needs zlo.type = surface_layer, which --terrain-method immersed does not use"),
    # the terrain tables hold one row per grid node
    (lambda s: s.terrain == "file_usgs" and (s.n[0] < 2 or s.n[1] < 2),
     "--terrain file_usgs needs at least 2 cells in x and y"),
    (lambda s: (s.case == "real") != (s.setup in REAL_SETUPS),
     "--case {case} does not match --setup {setup}; the real-data setups are wrfinput and metgrid"),
    (lambda s: s.real and s.init_given is not None,
     "--init is set by --setup {setup}; remove --init"),
    (lambda s: (s.lateral == "real") != s.real,
     "--lateral real belongs to the real-data setups (wrfinput, metgrid), and they take no other lateral type"),
    (lambda s: s.real and (s.terrain != "none" or s.terrain_method != "fitted" or s.stretched or s.buildings_on),
     "--setup {setup} takes its terrain and vertical levels from the NetCDF file; remove the terrain, "
     "--dz0 and --buildings options"),
    # ERF_DataStruct.H:1007-1010: WRFInput stops without a moisture model
    (lambda s: s.setup == "wrfinput" and not s.moist,
     "--setup wrfinput needs a moisture model; ERF stops with --moisture none"),
    (lambda s: s.real and s.forcing != "none",
     "--forcing {forcing} is an idealised driver; a real-data case is driven by its boundary file"),
    # ERF_DataStruct.H:1057-1063: the SAM schemes stop on a terrain-fitted or stretched mesh
    (lambda s: s.moisture.startswith("sam") and s.method == "fitted" and s.terrain != "none",
     "--moisture {moisture} does not run on a terrain-fitted mesh; use kessler, morrison, wsm6 or wdm6"),
    (lambda s: s.moisture.startswith("sam") and s.real,
     "--moisture {moisture} does not run on the terrain-fitted mesh of a real-data case"),
    # ERF prints that the land model "will not be applied" when no face is a surface layer
    (lambda s: s.lsm != "none" and not s.zlo_most,
     "--lsm {lsm} acts through zlo.type = surface_layer; add --surface most (not with --terrain-method immersed)"),
    # Source/LandSurfaceModel/SLM/ERF_SLM.H:211 asserts a moisture model
    (lambda s: s.lsm != "none" and not s.moist,
     "--lsm {lsm} needs a moisture model; ERF stops with --moisture none"),
    # ERF_CanopyTurbStruct.H / ERF.cpp:4229-4307
    (lambda s: s.canopy != "none" and s.method == "eb",
     "--canopy file is not covered with --terrain-method eb"),
    # Source/SourceTerms/ERF_MakeMomSources.cpp:119-121
    (lambda s: s.canopy != "none" and s.method == "immersed",
     "--canopy file cannot be used with --terrain-method immersed; ERF stops"),
    # ERF_SurfaceLayer.H:841: the reference height (first cell centre) must exceed the roughness length
    (lambda s: s.zlo_most and 0.5 * s.dz_first <= s.z0,
     "--z0 {z0} m is not below the first cell centre, {half_first} m; ERF stops. Lower --z0 or thicken the first cell"),
    # Source/LinearSolvers/ERF_PoissonWallDist.cpp: the wall-distance solve needs a no-slip or surface-layer face
    (lambda s: s.turbulence == "keqn" and s.surface == "slip" and s.method != "immersed",
     "--turbulence keqn needs a solid boundary for its wall distance; use --surface noslip or most"),
    # seen at ERF abce396eb: write_1D_profiles ends in a segmentation fault without a closure
    (lambda s: "profiles" in s.outputs and s.turbulence == "none",
     "--output profiles needs a turbulence closure; with --turbulence none ERF crashes at the first profile write"),
    # Source/IO/ERF_SampleData.H:2193: plane sampling over levels needs ref_ratio 2 2 2; this tool refines x and y
    (lambda s: "plane" in s.outputs and s.amr > 0,
     "--output plane cannot be combined with --amr; plane sampling over levels needs refinement in z as well"),
    (lambda s: s.lateral == "slipwall" and not s.at_rest,
     "--lateral slipwall closes the domain; a mean wind cannot pass through it. Use --wind 0 0 or another lateral type"),
    (lambda s: s.lateral == "inflow_outflow" and not s.wind[0] > 0.0,
     "--lateral inflow_outflow puts the inflow at xlo; give --wind with a positive first component"),
    (lambda s: s.setup == "density_current" and not s.at_rest,
     "--setup density_current starts at rest; the problem does not read a wind. Remove --wind"),
    # seen at ERF abce396eb: EB cases stopped at start-up when compressible ("Base state is too far out of HSE")
    # and diverged within 30 steps over a round hill; a ridge with the anelastic set ran 100 steps
    (lambda s: s.method == "eb" and (s.dycore != "anelastic" or not s.ridge),
     "--terrain-method eb is covered only for a ridge (fewer than 16 cells in y) with --dycore anelastic; "
     "other EB cases stopped or diverged in test runs"),
    # seen at ERF abce396eb over 100 steps: with --setup abl a ridge from a table diverged as an EB,
    # while the analytic ridge and the hill setup stayed steady
    (lambda s: s.method == "eb" and s.setup == "abl" and s.terrain in ("file", "file_usgs"),
     "--terrain-method eb with terrain from a table is covered for --setup hill only; with --setup abl it "
     "diverged in test runs. Use --terrain analytic, or --terrain-method fitted or immersed"),
    # the Schar ripples and the cos^4 hill must be resolved and fit in the domain
    (lambda s: s.terrain == "analytic" and s.hill == "schar" and 0.8 * s.hill_width < 8.0 * s.dx[0],
     "--hill schar: the ripple wavelength, 0.8 x --hill-width, must span at least 8 cells; raise --hill-width "
     "or --nx"),
    (lambda s: s.terrain == "analytic" and s.hill == "cos4" and s.ridge,
     "--hill cos4 is a round hill; give at least 16 cells in y (--ny), or use --hill woa for a ridge"),
    # the sea-surface table sets the roughness from the sea state
    (lambda s: s.zlo_most and s.surface_heat == "file" and s.roughness == "file",
     "--surface-heat file is a sea-surface temperature table over water, where the roughness follows the sea "
     "state; remove --roughness file"),
    (lambda s: s.terrain == "analytic" and s.hill == "cos4" and not s.ridge
     and 8.0 * s.hill_width > min(s.hi[0] - s.lo[0], s.hi[1] - s.lo[1]),
     "--hill cos4 ends at 4 x --hill-width from its centre; that must fit in the domain"),
    (lambda s: s.amr > 0 and s.method in ("eb",),
     "--amr with --terrain-method eb is not covered by this tool"),
    (lambda s: s.amr > 0 and s.refine == "box" and (s.n[0] % 4 or s.n[1] % 4),
     "--amr needs --nx and --ny divisible by 4 so the refined box lies on coarse cells"),
    (lambda s: s.damping == "sponge" and s.real,
     "--damping sponge is not covered for real-data cases; use rayleigh"),
    (lambda s: s.amr == 0 and (s.refine != "box" or s.refine_on is not None or s.refine_box is not None),
     "--refine, --refine-on and --refine-box need --amr 1 or 2"),
    # seen at ERF 329b9411c: the PBL kernels assert that every box spans the full column, and tagged boxes need not
    (lambda s: s.amr > 0 and s.refine == "field" and s.turbulence in PBL,
     "--refine field cannot be used with a PBL scheme: the scheme needs boxes over the full column, and tagged boxes "
     "need not reach the ground; use --refine box"),
    # seen at ERF 329b9411c: 3-D anelastic cases with field tagging stopped in the pressure solve or with a negative
    # theta within 20 steps, with and without a closure; the 2-D density current ran. After erf-model/ERF#5075
    # (7961a9f35) three of six such cases ran and three still stopped (the multigrid solve, an FFT boundary assertion)
    (lambda s: s.amr > 0 and s.refine == "field" and s.dycore == "anelastic" and s.n[1] > 1,
     "--refine field with --dycore anelastic is covered only for 2-D cases (--ny 1); half of the 3-D cases stopped "
     "in test runs. "
     "Use --dycore compressible or --refine box"),
    # Source/Prob/ERF_InitCustomPert_Bubble.H asserts that every box spans the full depth
    (lambda s: s.amr > 0 and s.refine == "field" and s.setup == "bubble",
     "--refine field cannot be used with --setup bubble: the Bubble problem needs boxes over the full depth; "
     "use --refine box"),
    (lambda s: s.amr > 0 and s.refine == "box" and (s.refine_on is not None or s.refine_box is not None),
     "--refine-on and --refine-box belong to --refine field"),
    # Source/Refinement/ERF_Tagging.cpp check_field: a moisture field needs a model that carries it
    (lambda s: any(f in ("qv", "qc", "qt") for f, _, _ in s.conditions) and not s.moist,
     "--refine-on: qv, qc and qt need a moisture model (--moisture)"),
    (lambda s: any(f == "qr" for f, _, _ in s.conditions) and s.moisture not in RAIN_MODELS,
     "--refine-on: qr needs a moisture model with rain; --moisture {moisture} has none"),
    (lambda s: any(f in ("qi", "qs", "qg", "max_reflectivity") for f, _, _ in s.conditions)
     and s.moisture not in ICE_MODELS,
     "--refine-on: qi, qs, qg and max_reflectivity need a moisture model with ice; --moisture {moisture} has none"),
    (lambda s: any(f == "terrain_blanking" for f, _, _ in s.conditions) and s.method != "immersed",
     "--refine-on: terrain_blanking needs --terrain-method immersed"),
    # ERF_Tagging.cpp:952-956: a box height above ground does not follow immersed or EB terrain
    (lambda s: s.refine_box is not None and len(s.refine_box) == 6 and s.refine_z == "above_ground"
     and s.method in ("immersed", "eb"),
     "--refine-z above_ground does not follow immersed or EB terrain; ERF stops. Use --refine-z physical"),
    # netcdf is the file type of the plotfile stream, not a stream of its own
    (lambda s: "netcdf" in s.outputs and "plt" not in s.outputs,
     "--output netcdf selects the file type of the plotfiles; add plt to --output"),
]


def check(s) -> List[str]:
    """Messages for every violated range or rule; empty when the case can be written."""
    errors = _ranges(s)
    if errors:
        return errors
    fields = dict(vars(s))
    fields["half_first"] = 0.5 * s.dz_first
    for test, message in RULES:
        if test(s):
            errors.append(message.format(**fields))
    return errors
