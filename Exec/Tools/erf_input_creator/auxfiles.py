"""Dummy auxiliary files: valid in format, placeholder in content.

Each writer returns an ``AuxFile`` whose text follows the layout the ERF
reader expects, sized to the grid of the case, so the run starts.  The values
are uniform placeholders that the user replaces.  ``describe`` holds what the
printed summary and FILES_TO_REPLACE.txt say about the file.
"""

from typing import List

import math

from .recipes.physics import GEO_TABLE_FILE, ROUGHNESS_FILE, SFC_FILE, SOUNDING_FILE
from .deckwriter import fmt_value as num
from .recipes.models import FOREST_FILE, RRTMGP_FILES, WINDFARM_LOC_FILE, WINDFARM_SPEC_FILE
from .recipes.physics import METEM_FILES, WRFBDY_FILE, WRFINPUT_FILE
from .recipes.terrain import BUILDINGS_FILE, TERRAIN_FILE, TERRAIN_USGS_FILE


class AuxFile:
    """A file the deck names.  ``text`` is None for a file the tool cannot generate."""

    def __init__(self, name: str, key: str, text, describe: List[str]):
        self.name = name
        self.key = key
        self.text = text
        self.describe = describe

    @property
    def generated(self) -> bool:
        return self.text is not None


def sounding(s) -> AuxFile:
    """erf.input_sounding_file; reader: Source/DataStructs/ERF_InputSoundingData.H."""
    zlo, zhi = s.lo[2], s.hi[2]
    nrow = min(max(s.n[2], 2), 50)
    # theta is written exactly as in the deck: ERF compares erf.most.pot_surf_temp with it for equality
    rows = ["1000.0 %s 0.0" % num(float(s.theta))]
    for k in range(nrow + 1):
        z = zhi if k == nrow else zlo + (zhi - zlo) * k / nrow      # the last row is the domain top exactly
        rows.append("%s %s 0.0 %s %s" % (num(float(z)), num(float(s.theta)), num(float(s.wind[0])), num(float(s.wind[1]))))
    return AuxFile(SOUNDING_FILE, "erf.input_sounding_file", "\n".join(rows) + "\n", [
        "Dummy content: neutral (theta = %g K), dry, uniform wind (%g, %g) m/s." % (s.theta, s.wind[0], s.wind[1]),
        "Format: no header, no comments, no blank lines.",
        "  line 1 : surface pressure [hPa], surface theta [K], surface qv [g/kg]",
        "  lines 2+: z [m], theta [K], qv [g/kg], u [m/s], v [m/s], with z strictly increasing",
        "A row at the domain bottom must repeat the surface theta and qv of line 1.",
        "The profile should reach the domain top (%g m); ERF extrapolates above the last row." % zhi,
        "erf.most.pot_surf_temp, when set, must equal the surface theta of line 1.",
    ])


def roughness(s) -> AuxFile:
    """erf.most.roughness_file_name; reader: Source/BoundaryConditions/ERF_SurfaceLayer.cpp."""
    nx, ny = s.n[0], s.n[1]
    rows = []
    for j in range(ny + 1):
        y = s.lo[1] + j * s.dx[1]
        for i in range(nx + 1):
            rows.append("%s %s %s" % (num(float(s.lo[0] + i * s.dx[0])), num(float(y)), num(float(s.z0))))
    return AuxFile(ROUGHNESS_FILE, "erf.most.roughness_file_name", "\n".join(rows) + "\n", [
        "Dummy content: z0 = %g m everywhere." % s.z0,
        "Format: no header; one row 'x y z0' [m] per node of the level-0 grid.",
        "  %d x %d nodes = %d rows; x runs fastest: x = %g + i*%g (i = 0..%d), y = %g + j*%g (j = 0..%d)"
        % (nx + 1, ny + 1, (nx + 1) * (ny + 1), s.lo[0], s.dx[0], nx, s.lo[1], s.dx[1], ny),
        "The x and y of each row must match the node positions; z0 must be positive and below erf.most.zref.",
        "Regenerate the file when amr.n_cell or the domain changes.",
    ])


def sfc_table(s) -> AuxFile:
    """erf.most.sfc_file; reader: read_cols in Source/BoundaryConditions/ERF_SurfaceLayer.cpp."""
    rows = ["day sst(K)", "0.0 %s" % num(float(s.theta)), "1.0 %s" % num(float(s.theta))]
    return AuxFile(SFC_FILE, "erf.most.sfc_file", "\n".join(rows) + "\n", [
        "Dummy content: constant sea-surface temperature %g K for one day." % s.theta,
        "ERF applies the table over water only; the deck sets erf.is_land = 0.",
        "Format: one header line, then one row per time.",
        "  2 columns (erf.most.use_sfc_sst)   : day, surface temperature [K]",
        "  5 columns (erf.most.use_sfc_fluxes): day, sst [K], H [W/m2], LE [W/m2], TAU [m2/s2]",
        "Time is in days; the first row is taken as the start of the run.",
    ])


def geo_table(s) -> AuxFile:
    """erf.abl_geo_wind_table; reader: Source/Initialization/ERF_InitGeowind.cpp."""
    rows = ["# z [m]   Ug [m/s]   Vg [m/s]",
            "%s %s %s" % (num(float(s.lo[2])), num(float(s.wind[0])), num(float(s.wind[1]))),
            "%s %s %s" % (num(float(s.hi[2])), num(float(s.wind[0])), num(float(s.wind[1])))]
    return AuxFile(GEO_TABLE_FILE, "erf.abl_geo_wind_table", "\n".join(rows) + "\n", [
        "Dummy content: geostrophic wind (%g, %g) m/s at every height." % (s.wind[0], s.wind[1]),
        "Format: rows 'z Ug Vg'; lines starting with # and blank lines are allowed.",
        "z must increase and span the domain from %g to %g m." % (s.lo[2], s.hi[2]),
    ])


def _nodes(s):
    # the last node is the domain face exactly, so the tables always cover the domain
    xs = [s.hi[0] if i == s.n[0] else s.lo[0] + i * s.dx[0] for i in range(s.n[0] + 1)]
    ys = [s.hi[1] if j == s.n[1] else s.lo[1] + j * s.dx[1] for j in range(s.n[1] + 1)]
    return xs, ys


def dummy_hill(s, x: float, y: float) -> float:
    """A cos^2 hill at the domain centre, 5 percent of the depth high, zero on the lateral faces.

    Zero at the faces keeps the table periodic.  With fewer than 16 cells in y it is a ridge across x.
    """
    lx, ly, lz = (s.hi[k] - s.lo[k] for k in range(3))
    rx = (x - 0.5 * (s.lo[0] + s.hi[0])) / (0.25 * lx)
    ry = (y - 0.5 * (s.lo[1] + s.hi[1])) / (0.25 * ly) if not s.ridge else 0.0
    r = math.sqrt(rx * rx + ry * ry)
    return 0.05 * lz * math.cos(0.5 * math.pi * r) ** 2 if r < 1.0 else 0.0


def dummy_building(s, x: float, y: float) -> float:
    """One block at the domain centre: a tenth of the depth high, an eighth of the extent wide."""
    lx, ly, lz = (s.hi[k] - s.lo[k] for k in range(3))
    inside_x = abs(x - 0.5 * (s.lo[0] + s.hi[0])) <= lx / 16.0 + 1e-9
    inside_y = abs(y - 0.5 * (s.lo[1] + s.hi[1])) <= ly / 16.0 + 1e-9 or s.ridge
    return 0.1 * lz if inside_x and inside_y else 0.0


def _table(s, height) -> str:
    """The ERF text table: nx, ny, x values, y values, heights with y fastest, one value per line."""
    xs, ys = _nodes(s)
    rows = ["%d" % len(xs), "%d" % len(ys)]
    rows += [num(float(x)) for x in xs] + [num(float(y)) for y in ys]
    rows += ["%.6f" % height(s, x, y) for x in xs for y in ys]
    return "\n".join(rows) + "\n"


_TABLE_FORMAT = [
    "Format: one number per line, no comments.",
    "  line 1: nx (number of x values); line 2: ny",
    "  next nx lines: x [m], strictly increasing; next ny lines: y [m], strictly increasing",
    "  next nx*ny lines: height [m], y fastest: value n belongs to (x[n // ny], y[n % ny])",
    "x and y are metres in the ERF domain, not longitude and latitude. ERF interpolates bilinearly,",
    "so the table need not match the grid.",
]


def terrain_table(s) -> AuxFile:
    """erf.terrain_file_name; reader: Source/Utils/ERF_TerrainFileChecks.cpp."""
    return AuxFile(TERRAIN_FILE, "erf.terrain_file_name", _table(s, dummy_hill), [
        "Dummy content: one smooth hill, %g m high, at the domain centre; zero at the lateral faces."
        % (0.05 * (s.hi[2] - s.lo[2])),
    ] + _TABLE_FORMAT + [
        "The table must cover the domain: each end may stop at most one table spacing inside it.",
        "In a periodic direction the heights on the two faces must be equal.",
    ])


def terrain_usgs(s) -> AuxFile:
    """erf.terrain_file_name_USGS; reader: Source/Utils/ERF_TerrainFileChecks.cpp."""
    xs, ys = _nodes(s)
    rows = ["0.0 0.0", "%d %d" % (len(xs), len(ys))]
    rows += ["%s %s %.6f" % (num(float(x)), num(float(y)), dummy_hill(s, x, y)) for y in ys for x in xs]
    return AuxFile(TERRAIN_USGS_FILE, "erf.terrain_file_name_USGS", "\n".join(rows) + "\n", [
        "Dummy content: one smooth hill, %g m high, at the domain centre; origin longitude and latitude 0, 0."
        % (0.05 * (s.hi[2] - s.lo[2])),
        "Format:",
        "  line 1: longitude and latitude [deg] of the table origin",
        "  line 2: nx ny (both at least 2)",
        "  next nx*ny lines: x [m], y [m], height [m], with x fastest; every row of one y repeats that y",
        "x increases along each row and y from row to row. The table must cover the domain.",
    ])


def buildings_table(s) -> AuxFile:
    """erf.buildings_file_name; same reader as the terrain table, without the coverage check."""
    return AuxFile(BUILDINGS_FILE, "erf.buildings_file_name", _table(s, dummy_building), [
        "Dummy content: one block, %g m high, at the domain centre." % (0.1 * (s.hi[2] - s.lo[2])),
    ] + _TABLE_FORMAT + [
        "Heights are building heights above the ground; zero means no building.",
        "The table may cover only part of the domain; outside it there are no buildings.",
    ])


def forest(s) -> AuxFile:
    """erf.forest_file; reader: Source/SourceTerms/ERF_ForestDrag.cpp."""
    lx, ly, lz = (s.hi[k] - s.lo[k] for k in range(3))
    xc, yc = 0.5 * (s.lo[0] + s.hi[0]), 0.5 * (s.lo[1] + s.hi[1])
    height, diameter = min(20.0, 0.1 * lz), 0.25 * min(lx, ly)
    row = "1 %s %s %s %s 0.2 6.0 0.8" % tuple(num(float(v)) for v in (xc, yc, height, diameter))
    return AuxFile(FOREST_FILE, "erf.forest_file", row + "\n", [
        "Dummy content: one patch of trees, %g m tall and %g m across, at the domain centre." % (height, diameter),
        "Format: no header; one row of 8 numbers per patch:",
        "  type (1 or 2), xc [m], yc [m], tree height [m], patch diameter [m], drag coefficient cd,",
        "  leaf area index lai, laimax (fraction of the height where the leaf area density peaks)",
        "A row with fewer than 8 numbers ends the read without a message.",
    ])


def windfarm_loc(s) -> AuxFile:
    """erf.windfarm_loc_table with windfarm_loc_type = x_y; reader: ERF_WindFarm.cpp."""
    xc, yc = 0.5 * (s.lo[0] + s.hi[0]), 0.5 * (s.lo[1] + s.hi[1])
    return AuxFile(WINDFARM_LOC_FILE, "erf.windfarm_loc_table", "%s %s\n" % (num(float(xc)), num(float(yc))), [
        "Dummy content: one turbine at the domain centre.",
        "Format: no header; one row 'x y' [m] per turbine, in the ERF domain.",
    ])


def windfarm_spec(s) -> AuxFile:
    """erf.windfarm_spec_table; reader: ERF_WindFarm.cpp (read_tables)."""
    lz = s.hi[2] - s.lo[2]
    hub = min(100.0, 0.25 * lz)
    rows = ["4", "%.1f %.1f 0.130 2.0" % (hub, hub), "3.0 0.80 0.05", "8.0 0.80 0.90", "12.0 0.60 2.00", "25.0 0.10 2.00"]
    return AuxFile(WINDFARM_SPEC_FILE, "erf.windfarm_spec_table", "\n".join(rows) + "\n", [
        "Dummy content: a turbine with hub height %g m and rotor diameter %g m, and a 4-row curve." % (hub, hub),
        "Format:",
        "  line 1: number of rows of the curve (at least 2)",
        "  line 2: hub height [m], rotor diameter [m], standing thrust coefficient, nominal power",
        "  then one row per wind speed: wind speed [m/s], thrust coefficient (at most 1), power",
        "Power uses the unit of the turbine table the values are taken from; the dummy numbers are not a real turbine.",
        "The rotor radius must not exceed the hub height.",
    ])


def required_files(s) -> List[AuxFile]:
    """Files the tool cannot generate: binary or model data the user must supply."""
    files = []
    if s.setup == "wrfinput":
        files.append(AuxFile(WRFINPUT_FILE, "erf.nc_init_file_0", None, [
            "A WRF wrfinput file (NetCDF) from real.exe or ideal.exe.",
            "Its west_east and south_north dimensions and DX, DY give amr.n_cell and the domain extent;",
            "bottom_top must be at least the z cell count. The run starts at the file's SIMULATION_START_DATE.",
            "A small example is Exec/RegTests/NoahMP_Ideal/wrfinput_ideal.cdl (ncgen -o wrfinput_d01 <file>).",
        ]))
        files.append(AuxFile(WRFBDY_FILE, "erf.nc_bdy_file", None, [
            "The WRF wrfbdy file (NetCDF) of the same domain: at least two evenly spaced times, and the",
            "_BXS/_BXE/_BYS/_BYE values and tendencies of U, V, T, QVAPOR, PH, MU. Its bdy_width must be",
            "at least erf.real_width, and its times must cover the run.",
        ]))
    if s.setup == "metgrid":
        for name in METEM_FILES:
            files.append(AuxFile(name, "erf.nc_init_file_0", None, [
                "A WPS met_em file (NetCDF). Give two or more in time order with a constant interval and put",
                "their names in erf.nc_init_file_0. DX, DY and the grid dimensions must match the deck.",
                "Examples: Exec/RegTests/MetGrid/met_em.d01.2022-06-18_00_00_00.nc and _06_00_00.nc.",
            ]))
    if s.radiation == "rrtmgp":
        for key, name in RRTMGP_FILES:
            files.append(AuxFile(name, key, None, [
                "RRTMGP coefficient file (NetCDF) from the rte-rrtmgp data set, in erf.rrtmgp_file_path.",
            ]))
    if s.lsm == "noahmp":
        files.append(AuxFile("namelist.erf", "(run directory)", None, [
            "Noah-MP namelist, group &NOAHLSM_OFFLINE, naming the land setup file of each level",
            "(ERF_SETUP_FILE_01..03). Template: Exec/RegTests/NoahMP_Ideal/namelist.erf.",
        ]))
        files.append(AuxFile("NoahmpTable.TBL", "(run directory)", None, [
            "Noah-MP parameter table. Copy Exec/RegTests/NoahMP_Ideal/NoahmpTable.TBL.",
        ]))
    return files


def aux_files(s) -> List[AuxFile]:
    """The files the deck of this case names."""
    files = []
    if s.init == "sounding":
        files.append(sounding(s))
    if s.zlo_most and s.roughness == "file":
        files.append(roughness(s))
    if s.zlo_most and s.surface_heat == "file":
        files.append(sfc_table(s))
    if s.forcing == "geo_table":
        files.append(geo_table(s))
    if s.terrain == "file":
        files.append(terrain_table(s))
    if s.terrain == "file_usgs":
        files.append(terrain_usgs(s))
    if s.buildings_on:
        files.append(buildings_table(s))
    if s.canopy == "file":
        files.append(forest(s))
    if s.windfarm != "none":
        files += [windfarm_loc(s), windfarm_spec(s)]
    return files + required_files(s)
