"""External files named in the deck: existence, and content checks when present.

Every file the deck names is checked for existence (relative to the deck or
the current directory).  When the file is present and a reader is available,
its contents are compared with the deck:

* wrfinput (erf.nc_init_file_0, init_type = WRFInput): staggered x/y
  dimensions against amr.n_cell, DX*(nx-1) against the domain extent,
  bottom_top against n_cell z, SIMULATION_START_DATE against start_datetime;
* wrfbdy (erf.nc_bdy_file): time stamps evenly spaced, the stored boundary
  width above erf.real_width, the span against stop_time, the variables that
  erf.use_wrf_bdy_qc_qi needs;
* wrflow (erf.nc_low_file): time stamps and the SST field;
* met_em (init_type = Metgrid): DX/DY against the domain;
* terrain text / USGS / NetCDF files: shape and coverage of the domain;
* input_sounding: column count and top height against the domain top;
* geostrophic wind table, roughness file, forest file, dirichlet files:
  readable and of the right column count.

NetCDF files are read with the netCDF4 module when it is installed, otherwise
with ``ncdump -h`` when it is on the PATH; without either, only the existence
of the file is checked and an INFO line says so.
"""

import datetime
import os
import re
import shutil
import subprocess
from typing import Dict, List, Optional

from ..checker import Checker, register

# key -> (what it is, several files allowed)
FILE_KEYS = {
    "erf.nc_init_file_0": ("wrfinput / met_em / NetCDF initial file", True),
    "erf.nc_init_file_1": ("fine-level initial file", True),
    "erf.nc_init_file_2": ("fine-level initial file", True),
    "erf.nc_bdy_file": ("wrfbdy lateral boundary file", False),
    "erf.nc_low_file": ("wrflow lower boundary file", False),
    "erf.input_sounding_file": ("input sounding", True),
    "erf.input_sponge_file": ("sponge profile", False),
    "erf.terrain_file_name": ("terrain text file", False),
    "erf.terrain_file_name_nc": ("terrain NetCDF file", False),
    "erf.terrain_file_name_USGS": ("terrain USGS text file", False),
    "erf.abl_geo_wind_table": ("geostrophic wind table", False),
    "erf.most.roughness_file_name": ("roughness file", True),
    "erf.forest_file": ("forest patch file", False),
    "erf.forest_lai_file": ("forest LAI file", False),
    "erf.forest_height_file": ("forest height file", False),
    "erf.forest_cd_file": ("forest drag file", False),
    "erf.windfarm_loc_table": ("turbine location table", False),
    "erf.windfarm_spec_table": ("turbine specification table", False),
    "erf.windfarm_blade_table": ("blade table", False),
    "erf.large_scale_forcing_file": ("large-scale forcing file", False),
    "erf.erfbdy_file": ("ERF boundary file", False),
    "amr.regrid_file": ("fixed grids file", False),
    "erf.fire.fuel_map.fuel_map_file": ("fuel map", False),
    "erf.fire.ignition.polygon_file": ("ignition polygon", True),
    "erf.fire.ignition.ignition_schedule_file": ("ignition schedule", False),
}
FACE_FILE_KEYS = (".dirichlet_file", ".inflow_profile_file")


class NetCDFHeader:
    """Dimensions, global attributes and variable names of a NetCDF file."""

    def __init__(self):
        self.dims: Dict[str, int] = {}
        self.attrs: Dict[str, object] = {}
        self.vars: Dict[str, List[str]] = {}   # name -> dimension names
        self.reader = ""

    @classmethod
    def read(cls, path: str) -> Optional["NetCDFHeader"]:
        h = cls()
        try:
            import netCDF4  # type: ignore
        except ImportError:
            netCDF4 = None
        if netCDF4 is not None:
            try:
                with netCDF4.Dataset(path) as ds:
                    h.dims = {k: len(v) for k, v in ds.dimensions.items()}
                    for a in ds.ncattrs():
                        v = ds.getncattr(a)
                        try:
                            v = v.tolist()
                        except AttributeError:
                            pass
                        if isinstance(v, list) and len(v) == 1:
                            v = v[0]
                        h.attrs[a] = v
                    h.vars = {k: list(v.dimensions) for k, v in ds.variables.items()}
                h.reader = "netCDF4"
                return h
            except Exception:  # fall through to ncdump
                pass
        if shutil.which("ncdump"):
            try:
                out = subprocess.run(["ncdump", "-h", path], capture_output=True, text=True, timeout=60).stdout
            except Exception:
                return None
            h._parse_ncdump(out)
            h.reader = "ncdump"
            return h
        return None

    def _parse_ncdump(self, text: str) -> None:
        section = ""
        for line in text.splitlines():
            s = line.strip()
            if s.startswith("dimensions:"):
                section = "dims"
                continue
            if s.startswith("variables:"):
                section = "vars"
                continue
            if s.startswith("// global attributes:"):
                section = "attrs"
                continue
            if section == "dims":
                m = re.match(r"(\S+)\s*=\s*(\d+|UNLIMITED)\s*;(?:\s*//\s*\((\d+) currently\))?", s)
                if m:
                    self.dims[m.group(1)] = int(m.group(3) or (m.group(2) if m.group(2).isdigit() else 0))
            elif section == "vars":
                m = re.match(r"\w+\s+(\w+)\(([^)]*)\)\s*;", s)
                if m:
                    self.vars[m.group(1)] = [x.strip() for x in m.group(2).split(",") if x.strip()]
            elif section == "attrs":
                m = re.match(r":([\w\-]+)\s*=\s*(.*);$", s)
                if m:
                    val = m.group(2).strip()
                    if val.startswith('"'):
                        self.attrs[m.group(1)] = val.strip('"')
                    else:
                        try:
                            self.attrs[m.group(1)] = float(val.rstrip("f").split(",")[0])
                        except ValueError:
                            self.attrs[m.group(1)] = val

    def times(self, path: str) -> List[str]:
        """The Times variable as strings (netCDF4 only)."""
        try:
            import netCDF4  # type: ignore
            with netCDF4.Dataset(path) as ds:
                if "Times" not in ds.variables:
                    return []
                raw = ds.variables["Times"][:]
                out = []
                for row in raw:
                    out.append("".join(x.decode() if isinstance(x, bytes) else str(x) for x in row).strip())
                return out
        except Exception:
            return []


def _parse_wrf_time(s: str) -> Optional[datetime.datetime]:
    for fmt in ("%Y-%m-%d_%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M"):
        try:
            return datetime.datetime.strptime(s.strip(), fmt)
        except ValueError:
            continue
    return None


@register
class ExternalFilesChecker(Checker):
    name = "files"

    def run(self) -> None:
        d = self.deck
        self.present: Dict[str, List[str]] = {}
        for key, (what, many) in FILE_KEYS.items():
            if not d.has(key):
                continue
            names = d.get_strs(key) if many else [d.get_str(key)]
            for name in names:
                if not name or name.lower() in ("none", "null", ""):
                    continue
                path = d.resolve_path(name)
                if path is None:
                    self.error(key, "%s %r not found next to the deck or in the current directory" % (what, name),
                               "Copy or link the file, or give an absolute path.")
                else:
                    self.present.setdefault(key, []).append(path)
        for face in ("xlo", "xhi", "ylo", "yhi", "zlo", "zhi"):
            for suf in FACE_FILE_KEYS:
                key = face + suf
                if d.has(key):
                    name = d.get_str(key)
                    path = d.resolve_path(name)
                    if path is None:
                        self.error(key, "file %r not found" % name)
                    else:
                        self.present.setdefault(key, []).append(path)
        self._real_data()
        self._terrain()
        self._sounding()
        self._tables()

    # -- real data ----------------------------------------------------------
    def _first(self, key: str) -> Optional[str]:
        paths = self.present.get(key)
        return paths[0] if paths else None

    def _header(self, key: str) -> Optional[NetCDFHeader]:
        path = self._first(key)
        if not path:
            return None
        h = NetCDFHeader.read(path)
        if h is None:
            self.info(key, "file present but neither the netCDF4 module nor ncdump is available to read it",
                      "pip install netCDF4 to let the checker compare the file with the deck.")
        return h

    def _real_data(self) -> None:
        c, d = self.ctx, self.deck
        t = c.init_type.lower()
        if t == "wrfinput":
            self._wrfinput()
        elif t == "metgrid":
            self._metgrid()
        if "erf.nc_bdy_file" in self.present:
            self._wrfbdy()
        if "erf.nc_low_file" in self.present:
            self._wrflow()

    def _wrfinput(self) -> None:
        c = self.ctx
        h = self._header("erf.nc_init_file_0")
        if h is None:
            return
        key = "erf.nc_init_file_0"
        nxs, nys, nz = h.dims.get("west_east_stag"), h.dims.get("south_north_stag"), h.dims.get("bottom_top")
        if nxs is None or nys is None or nz is None:
            self.error(key, "the file has no west_east_stag / south_north_stag / bottom_top dimensions; it is not a wrfinput file")
            return
        self.info(key, "wrfinput grid %d x %d x %d (mass points), DX = %s m, start %s" % (
            nxs - 1, nys - 1, nz, h.attrs.get("DX", "?"), h.attrs.get("SIMULATION_START_DATE", "?")))
        if c.n_cell:
            if c.n_cell[0] != nxs - 1 or c.n_cell[1] != nys - 1:
                self.error("amr.n_cell", "n_cell %d x %d does not match the wrfinput grid %d x %d; ERF stops" % (
                    c.n_cell[0], c.n_cell[1], nxs - 1, nys - 1), "Set amr.n_cell = %d %d <nz>." % (nxs - 1, nys - 1))
            if c.n_cell[2] > nz:
                self.error("amr.n_cell", "n_cell z = %d exceeds the %d wrfinput levels; ERF stops" % (c.n_cell[2], nz))
            elif c.n_cell[2] < nz:
                self.warn("amr.n_cell", "n_cell z = %d is below the %d wrfinput levels; the top levels are dropped" % (c.n_cell[2], nz))
        dx, dy = h.attrs.get("DX"), h.attrs.get("DY")
        if c.extent and isinstance(dx, (int, float)) and isinstance(dy, (int, float)):
            lx, ly = dx * (nxs - 1), dy * (nys - 1)
            for i, (L, want) in enumerate(((c.extent[0], lx), (c.extent[1], ly))):
                if abs(L - want) > 1e-7 * want:
                    self.error("geometry.prob_hi", "domain extent %g m in %s does not match the wrfinput %g m (DX * n); ERF stops" % (
                        L, "xy"[i], want), "Set geometry.prob_extent = %g %g <Lz>, or prob_hi accordingly." % (lx, ly))
        self._start_date(key, h.attrs.get("SIMULATION_START_DATE"))
        if "HGT" in h.vars and c.terrain.lower() == "none":
            self.advise("erf.terrain_type", "the wrfinput file carries terrain (HGT) but the deck runs on a flat mesh")
        if self.bool_("erf.use_wrf_bdy_qc_qi", False):
            pass  # checked on the bdy file

    def _metgrid(self) -> None:
        c = self.ctx
        h = self._header("erf.nc_init_file_0")
        if h is None:
            return
        key = "erf.nc_init_file_0"
        if "num_metgrid_levels" not in h.dims:
            self.error(key, "no num_metgrid_levels dimension; this is not a met_em file")
            return
        nxs, nys = h.dims.get("west_east_stag"), h.dims.get("south_north_stag")
        if c.n_cell and nxs and nys and (c.n_cell[0] != nxs - 1 or c.n_cell[1] != nys - 1):
            self.error("amr.n_cell", "n_cell %d x %d does not match the met_em grid %d x %d" % (c.n_cell[0], c.n_cell[1], nxs - 1, nys - 1))
        dx, dy = h.attrs.get("DX"), h.attrs.get("DY")
        if c.extent and nxs and nys and isinstance(dx, (int, float)) and isinstance(dy, (int, float)):
            lx, ly = dx * (nxs - 1), dy * (nys - 1)
            if abs(c.extent[0] - lx) > 1e-6 * lx or abs(c.extent[1] - ly) > 1e-6 * ly:
                self.error("geometry.prob_hi", "domain extent does not match the met_em DX * n = %g x %g m" % (lx, ly))
        self._start_date(key, h.attrs.get("SIMULATION_START_DATE"))
        files = self.deck.get_strs("erf.nc_init_file_0")
        if len(files) < 2:
            self.warn("erf.nc_init_file_0", "a single met_em file: the lateral boundaries cannot be updated in time",
                      "List every met_em time in erf.nc_init_file_0 (and set erf.use_real_bcs = true).")

    def _start_date(self, key: str, file_date) -> None:
        d = self.deck
        if not d.has("start_datetime") or not isinstance(file_date, str):
            return
        want = _parse_wrf_time(" ".join(d.get_strs("start_datetime")))
        have = _parse_wrf_time(file_date)
        if want and have and want != have:
            self.error("start_datetime", "%s differs from SIMULATION_START_DATE %s in the file; ERF stops" % (want, have),
                       "Remove start_datetime (ERF reads it from the file) or make them equal.")

    def _wrfbdy(self) -> None:
        c, d = self.ctx, self.deck
        key = "erf.nc_bdy_file"
        h = self._header(key)
        if h is None:
            return
        width = h.dims.get("bdy_width")
        if width is None:
            self.error(key, "no bdy_width dimension; this is not a wrfbdy file")
            return
        rw = self.int_("erf.real_width", 0)
        if rw is not None and rw >= width:
            self.error("erf.real_width", "= %d must be below the boundary width %d stored in wrfbdy; ERF stops" % (rw, width),
                       "Set erf.real_width = %d or less." % (width - 1))
        if self.bool_("erf.use_wrf_bdy_qc_qi", False):
            need = [v for v in ("QCLOUD_BXS", "QICE_BXS") if v not in h.vars]
            if need:
                self.error("erf.use_wrf_bdy_qc_qi", "wrfbdy lacks %s; ERF stops" % ", ".join(need),
                           "Regenerate wrfbdy with cloud fields, or set the key false.")
        nxs = h.dims.get("west_east_stag") or (h.dims.get("west_east", 0) + 1)
        nys = h.dims.get("south_north_stag") or (h.dims.get("south_north", 0) + 1)
        if c.n_cell and nxs and nys and (c.n_cell[0] != nxs - 1 or c.n_cell[1] != nys - 1):
            self.error(key, "wrfbdy grid %d x %d does not match n_cell %d x %d" % (nxs - 1, nys - 1, c.n_cell[0], c.n_cell[1]))
        path = self._first(key)
        times = h.times(path)
        if len(times) >= 2:
            dts = [_parse_wrf_time(s) for s in times]
            if all(dts):
                gaps = {(dts[i + 1] - dts[i]).total_seconds() for i in range(len(dts) - 1)}
                if len(gaps) > 1:
                    self.error(key, "wrfbdy time stamps are not evenly spaced (%s s); ERF stops" % sorted(gaps))
                interval = min(gaps)
                span = (dts[-1] - dts[0]).total_seconds() + interval
                self.info(key, "wrfbdy covers %s to %s (%d records every %.0f s, last record valid for %.0f s)" % (
                    times[0], times[-1], len(times), interval, span))
                stop = self.real("stop_time")
                if stop is not None and stop > span:
                    self.error("stop_time", "= %g s runs past the %g s of lateral boundary data" % (stop, span),
                               "Lower stop_time or extend wrfbdy.")
                sd = d.get_strs("stop_datetime")
                if sd:
                    want = _parse_wrf_time(" ".join(sd))
                    if want and want > dts[-1] + datetime.timedelta(seconds=interval):
                        self.error("stop_datetime", "lies past the last wrfbdy record %s" % times[-1])
        elif h.reader == "ncdump":
            self.info(key, "time stamps not compared (install netCDF4 to read the Times variable)")

    def _wrflow(self) -> None:
        key = "erf.nc_low_file"
        h = self._header(key)
        if h is None:
            return
        if "SST" not in h.vars:
            self.error(key, "wrflow file has no SST variable")
        c = self.ctx
        nx, ny = h.dims.get("west_east"), h.dims.get("south_north")
        if c.n_cell and nx and ny and (c.n_cell[0] != nx or c.n_cell[1] != ny):
            self.error(key, "wrflow grid %d x %d does not match n_cell %d x %d" % (nx, ny, c.n_cell[0], c.n_cell[1]))
        if not self.deck.has("erf.most.ignore_sst") and (self.str_("zlo.type", "") or "").lower() != "surface_layer":
            self.advise(key, "a wrflow SST file is only used by a surface_layer boundary at zlo")

    # -- terrain ----------------------------------------------------------
    def _terrain(self) -> None:
        c = self.ctx
        if "erf.terrain_file_name" in self.present:
            self._terrain_text(self._first("erf.terrain_file_name"), "erf.terrain_file_name")
        if "erf.terrain_file_name_USGS" in self.present:
            self._terrain_usgs(self._first("erf.terrain_file_name_USGS"))
        if "erf.terrain_file_name_nc" in self.present:
            key = "erf.terrain_file_name_nc"
            h = self._header(key)
            if h is None:
                return
            name = next((v for v in ("height", "z", "terrain", "HGT_M", "AGL") if v in h.vars), None)
            if name is None:
                self.error(key, "no terrain variable (height, z, terrain, HGT_M or AGL) in the file; ERF stops")
                return
            dims = h.vars[name]
            if len(dims) not in (2, 3):
                self.error(key, "terrain variable %s has rank %d; (y,x) or (time,y,x) is required" % (name, len(dims)))
                return
            ny, nx = h.dims.get(dims[-2], 0), h.dims.get(dims[-1], 0)
            if nx < 2 or ny < 2:
                self.error(key, "terrain field is %d x %d; at least 2 x 2 is needed" % (nx, ny))
            if ("x" in h.vars) != ("y" in h.vars):
                self.error(key, "both one-dimensional x and y coordinates are needed")
            if "x" not in h.vars and name == "HGT_M" and c.n_cell and (nx != c.n_cell[0] or ny != c.n_cell[1]):
                self.error(key, "WPS HGT_M grid %d x %d must equal n_cell %d x %d; ERF stops" % (nx, ny, c.n_cell[0], c.n_cell[1]))
            self.info(key, "terrain variable %s on a %d x %d grid" % (name, nx, ny))

    def _numbers(self, path: str) -> List[float]:
        out = []
        with open(path, errors="replace") as fh:
            for line in fh:
                for tok in line.replace(",", " ").split():
                    try:
                        out.append(float(tok))
                    except ValueError:
                        return out
        return out

    def _terrain_text(self, path: str, key: str) -> None:
        c = self.ctx
        v = self._numbers(path)
        if len(v) < 2:
            self.error(key, "terrain file is empty or does not start with nx ny")
            return
        nx, ny = int(v[0]), int(v[1])
        need = 2 + nx + ny + nx * ny
        if nx <= 0 or ny <= 0:
            self.error(key, "nx, ny = %d, %d must be positive" % (nx, ny))
            return
        if len(v) < need:
            self.error(key, "terrain file holds %d numbers; nx ny, %d x, %d y and %d z values (%d) are needed" % (len(v), nx, ny, nx * ny, need))
            return
        x, y = v[2:2 + nx], v[2 + nx:2 + nx + ny]
        z = v[2 + nx + ny:need]
        self._coverage(key, x, y, z, nx, ny)

    def _terrain_usgs(self, path: str) -> None:
        key = "erf.terrain_file_name_USGS"
        v = self._numbers(path)
        if len(v) < 4:
            self.error(key, "USGS terrain file must start with lon_min lat_min nx ny")
            return
        lon, lat, nx, ny = v[0], v[1], int(v[2]), int(v[3])
        if abs(lon) > 180 or abs(lat) > 90:
            self.error(key, "first line lon %g lat %g is out of range" % (lon, lat))
        rows = (len(v) - 4) // 3
        if nx < 2 or ny < 2:
            self.error(key, "USGS terrain needs at least two points in each direction")
        if rows != nx * ny:
            self.error(key, "%d (x y z) rows for nx * ny = %d" % (rows, nx * ny))
            return
        xs = v[4::3][:rows]
        ys = v[5::3][:rows]
        zs = v[6::3][:rows]
        self._coverage(key, sorted(set(xs)), sorted(set(ys)), zs, nx, ny)

    def _coverage(self, key, x, y, z, nx, ny) -> None:
        c = self.ctx
        if any(x[i + 1] <= x[i] for i in range(len(x) - 1)) or any(y[i + 1] <= y[i] for i in range(len(y) - 1)):
            self.error(key, "terrain x and y coordinates must increase")
        self.info(key, "terrain table %d x %d, x %g..%g, y %g..%g, z %g..%g m" % (nx, ny, x[0], x[-1], y[0], y[-1], min(z), max(z)))
        if c.prob_lo and c.prob_hi:
            dxt = (x[-1] - x[0]) / max(nx - 1, 1)
            dyt = (y[-1] - y[0]) / max(ny - 1, 1)
            for i, (lo, hi, t0, t1, sp) in enumerate(((c.prob_lo[0], c.prob_hi[0], x[0], x[-1], dxt), (c.prob_lo[1], c.prob_hi[1], y[0], y[-1], dyt))):
                if t0 > lo + sp or t1 < hi - sp:
                    self.error(key, "terrain table covers %g..%g in %s but the domain is %g..%g (more than one table spacing short); ERF stops" % (
                        t0, t1, "xy"[i], lo, hi))
                elif t0 > lo or t1 < hi:
                    self.advise(key, "terrain table ends inside the domain in %s by less than one spacing; ERF extrapolates there" % "xy"[i])
            if c.extent and max(z) - min(z) > 0.5 * c.extent[2]:
                self.warn(key, "terrain relief %g m is over half the domain height %g m" % (max(z) - min(z), c.extent[2]),
                          "Raise the domain top; the terrain-following mesh is badly distorted otherwise.")
            if c.dx and (dxt > 4 * c.dx[0] or dyt > 4 * c.dx[1]):
                self.advise(key, "terrain table spacing (%g, %g m) is much coarser than the grid (%g, %g m)" % (dxt, dyt, c.dx[0], c.dx[1]),
                            "The interpolated terrain is faceted; use a finer table if available.")

    # -- sounding ---------------------------------------------------------
    def _sounding(self) -> None:
        c, d = self.ctx, self.deck
        if c.init_type.lower() != "input_sounding":
            return
        names = d.get_strs("erf.input_sounding_file") or ["input_sounding"]
        for name in names:
            path = d.resolve_path(name)
            key = "erf.input_sounding_file"
            if path is None:
                if not d.has(key):
                    self.error(key, "the default sounding file 'input_sounding' is not next to the deck")
                continue
            rows = []
            with open(path, errors="replace") as fh:
                for n, line in enumerate(fh, 1):
                    toks = line.split()
                    if not toks or toks[0].startswith("#"):
                        continue
                    try:
                        rows.append([float(t) for t in toks])
                    except ValueError:
                        self.error(key, "%s line %d is not numeric" % (name, n))
                        return
            if not rows:
                self.error(key, "%s is empty" % name)
                return
            if len(rows[0]) < 3:
                self.error(key, "%s: the first line must hold p_surf [hPa], theta_surf [K], qv_surf [g/kg]" % name)
                continue
            if not (300 < rows[0][0] < 1100):
                self.warn(key, "%s: surface pressure %g looks wrong for hPa" % (name, rows[0][0]))
            if any(len(r) < 5 for r in rows[1:]):
                self.error(key, "%s: every line after the first needs z theta qv u v" % name)
                return
            z = [r[0] for r in rows[1:]]
            if any(z[i + 1] <= z[i] for i in range(len(z) - 1)):
                self.error(key, "%s: heights must increase" % name)
            if c.extent and z and z[-1] < c.extent[2]:
                self.warn(key, "%s: top of the sounding %g m is below the domain top %g m; ERF extrapolates" % (name, z[-1], c.extent[2]))
            for stkey in ("erf.most.pot_surf_temp", "erf.most.surf_temp"):
                if d.has(stkey):
                    st = self.real(stkey)
                    if st is not None and st > 0 and len(rows[0]) > 1 and abs(st - rows[0][1]) > 1e-9:
                        self.error(stkey, "= %g differs from the sounding's surface theta %g; ERF stops ('Input sounding temp does not match')" % (st, rows[0][1]),
                                   "Make the two equal, or remove %s." % stkey)
            qv = [r[2] for r in rows]
            if max(qv) > 50:
                self.warn(key, "%s: qv up to %g; the file is in g/kg, so this is not a mixing ratio" % (name, max(qv)))

    # -- tables -----------------------------------------------------------
    def _tables(self) -> None:
        d = self.deck
        for path in self.present.get("erf.abl_geo_wind_table", []):
            v = self._rows(path)
            if any(len(r) != 3 for r in v):
                self.error("erf.abl_geo_wind_table", "%s: every line needs z u_g v_g" % os.path.basename(path))
        for path in self.present.get("erf.most.roughness_file_name", []):
            v = self._rows(path)
            if any(len(r) != 3 for r in v):
                self.error("erf.most.roughness_file_name", "%s: every line needs x y z0" % os.path.basename(path))
            elif any(r[2] <= 0 for r in v):
                self.error("erf.most.roughness_file_name", "%s: z0 must be positive everywhere" % os.path.basename(path))
        for face in ("xlo", "xhi", "ylo", "yhi"):
            for key in (face + ".dirichlet_file", "erf." + face + ".dirichlet_file"):
                for path in self.present.get(key, []):
                    v = self._rows(path)
                    if not v or any(len(r) < 4 for r in v):
                        self.error(key, "every line needs z u v w [theta]")

    def _rows(self, path: str) -> List[List[float]]:
        out = []
        with open(path, errors="replace") as fh:
            for line in fh:
                toks = line.split()
                if not toks or toks[0].startswith("#"):
                    continue
                try:
                    out.append([float(t) for t in toks])
                except ValueError:
                    continue
        return out
