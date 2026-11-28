"""Terrain checks: slope at each level's cell size, relief, water fraction, and the mismatch between
HRRR-like (3 km) terrain and the level-0 terrain along the lateral boundaries.

Elevation sources:

* ``terrarium`` (default): the public AWS Terrain Tiles (Mapzen / Tilezen, terrarium encoding,
  elevation = R*256 + G + B/256 - 32768 m), fetched over HTTPS and cached.  See the README for
  the data providers and the attribution they require.
* an ESRI ASCII grid file (``.asc``) in geographic coordinates (degrees), for example exported
  with ``gdal_translate -of AAIGrid`` from a DEM in EPSG:4326.
* ``none``: no terrain check.

The checks need numpy.  The slope is the largest central-difference gradient of the cell-mean
terrain on each level's grid; cell means come from k x k samples per cell.  The thresholds are
best guesses; see README.md.
"""

import http.client
import math
import os
import struct
import urllib.request
import zlib

from .plan import damping_depth

TILE_URL = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"
MAX_TILES = 400
MAX_SAMPLES = 4_000_000


class TerrainError(Exception):
    pass


# --------------------------------------------------------------------------------------------
# PNG decoding (8-bit greyscale, RGB or RGBA, not interlaced), enough for terrarium tiles.
# --------------------------------------------------------------------------------------------

def decode_png(data):
    """Return a numpy uint8 array (height, width, channels).  A damaged file raises TerrainError."""
    try:
        return _decode_png(data)
    except (zlib.error, struct.error, IndexError, ValueError) as e:
        raise TerrainError("damaged PNG data (%s)" % e)


def _decode_png(data):
    import numpy as np
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise TerrainError("not a PNG file")
    pos, idat, width = 8, [], None
    while pos < len(data):
        (length,) = struct.unpack(">I", data[pos:pos + 4])
        ctype = data[pos + 4:pos + 8]
        chunk = data[pos + 8:pos + 8 + length]
        pos += 12 + length
        if ctype == b"IHDR":
            width, height, depth, color, _, _, interlace = struct.unpack(">IIBBBBB", chunk)
            if depth != 8 or color not in (0, 2, 6) or interlace != 0:
                raise TerrainError("PNG layout not supported (depth %d, colour %d, interlace %d)"
                                   % (depth, color, interlace))
            nch = {0: 1, 2: 3, 6: 4}[color]
        elif ctype == b"IDAT":
            idat.append(chunk)
        elif ctype == b"IEND":
            break
    if width is None:
        raise TerrainError("PNG without IHDR")
    raw = zlib.decompress(b"".join(idat))
    stride = width * nch
    if len(raw) != height * (stride + 1):
        raise TerrainError("PNG data of %d bytes, %d expected" % (len(raw), height * (stride + 1)))
    out = bytearray(height * stride)
    prev = bytearray(stride)
    for r in range(height):
        f = raw[r * (stride + 1)]
        cur = bytearray(raw[r * (stride + 1) + 1:(r + 1) * (stride + 1)])
        if f == 1:
            for c in range(nch, stride):
                cur[c] = (cur[c] + cur[c - nch]) & 255
        elif f == 2:
            cur = bytearray((a + b) & 255 for a, b in zip(cur, prev))
        elif f == 3:
            for c in range(stride):
                left = cur[c - nch] if c >= nch else 0
                cur[c] = (cur[c] + ((left + prev[c]) >> 1)) & 255
        elif f == 4:
            for c in range(stride):
                if c >= nch:
                    a, cc = cur[c - nch], prev[c - nch]
                else:
                    a = cc = 0
                b = prev[c]
                pa, pb, pc = abs(b - cc), abs(a - cc), abs(a + b - 2 * cc)
                pred = a if (pa <= pb and pa <= pc) else (b if pb <= pc else cc)
                cur[c] = (cur[c] + pred) & 255
        elif f != 0:
            raise TerrainError("PNG filter %d unknown" % f)
        out[r * stride:(r + 1) * stride] = cur
        prev = cur
    return np.frombuffer(bytes(out), dtype=np.uint8).reshape(height, width, nch)


# --------------------------------------------------------------------------------------------
# Elevation sources.  Each has sample(lat, lon) -> elevation [m] on numpy arrays, and
# resolution(lat) -> its grid spacing [m].
# --------------------------------------------------------------------------------------------

class AsciiGridDEM:
    """ESRI ASCII grid in degrees: ncols, nrows, xllcorner|xllcenter, yllcorner|yllcenter,
    cellsize, optional NODATA_value, then rows from north to south."""

    def __init__(self, path):
        import numpy as np
        head = {}
        with open(path) as f:
            for _ in range(6):
                pos = f.tell()
                parts = f.readline().split()
                if len(parts) != 2 or parts[0][0].isdigit() or parts[0][0] == "-":
                    f.seek(pos)
                    break
                head[parts[0].lower()] = float(parts[1])
            data = np.loadtxt(f, dtype=float)
        for k in ("ncols", "nrows", "cellsize"):
            if k not in head:
                raise TerrainError("%s: header has no %s" % (path, k))
        self.nx, self.ny, self.cs = int(head["ncols"]), int(head["nrows"]), head["cellsize"]
        if data.shape != (self.ny, self.nx):
            raise TerrainError("%s: %d x %d values for a %d x %d header"
                               % (path, data.shape[0], data.shape[1], self.ny, self.nx))
        if "xllcenter" in head:
            self.lon0, self.lat0 = head["xllcenter"], head["yllcenter"]
        else:
            self.lon0, self.lat0 = head["xllcorner"] + self.cs / 2, head["yllcorner"] + self.cs / 2
        nod = head.get("nodata_value")
        self.z = data[::-1].copy()       # south to north
        if nod is not None:
            self.z[self.z == nod] = np.nan
        self.name = os.path.basename(path)

    def resolution(self, lat):
        return self.cs * 111195.0 * max(math.cos(math.radians(lat)), 0.1)

    def covers(self, lat, lon):
        lat1 = self.lat0 + (self.ny - 1) * self.cs
        lon1 = self.lon0 + (self.nx - 1) * self.cs
        return (lat.min() >= self.lat0 and lat.max() <= lat1 and
                lon.min() >= self.lon0 and lon.max() <= lon1)

    def sample(self, lat, lon):
        return _bilinear(self.z, (lon - self.lon0) / self.cs, (lat - self.lat0) / self.cs)


class TerrariumDEM:
    """AWS Terrain Tiles, terrarium PNG encoding, web-Mercator tiles of 256 x 256 pixels."""

    def __init__(self, cache_dir=None, fetch=None):
        self.cache = cache_dir or os.path.join(os.path.expanduser("~"), ".cache",
                                               "erf_domain_advisor", "terrarium")
        self.fetch = fetch or _http_get
        self.zoom = 10
        self.name = "AWS Terrain Tiles (terrarium)"
        self.mosaic = None
        self.despiked = 0
        self.lon_ref = 0.0

    def _unwrap(self, lon):
        """Longitudes within 180 deg of lon_ref, so a range across the antimeridian stays
        contiguous; tile columns past the last one wrap in prepare."""
        import numpy as np
        return self.lon_ref + (np.asarray(lon, dtype=float) - self.lon_ref + 180.0) % 360.0 - 180.0

    def resolution(self, lat):
        return 156543.03 * math.cos(math.radians(lat)) / 2 ** self.zoom

    def prepare(self, lat, lon, dx):
        """Choose the zoom (pixel at most dx / 2, zoom 15 at most) and load the tiles over the
        lat/lon range; drops zoom while more than MAX_TILES tiles would be needed."""
        import numpy as np
        latc = float(np.mean(lat))
        lr = np.radians(np.asarray(lon, dtype=float))
        self.lon_ref = float(np.degrees(np.arctan2(np.mean(np.sin(lr)), np.mean(np.cos(lr)))))
        lon = self._unwrap(lon)
        z = 0
        while z < 15 and 156543.03 * math.cos(math.radians(latc)) / 2 ** z > dx / 2:
            z += 1
        while True:
            x0, y1 = _tile_xy(lat.min(), lon.min(), z)
            x1, y0 = _tile_xy(lat.max(), lon.max(), z)
            tx0, tx1, ty0, ty1 = int(x0), int(x1), int(y0), int(y1)
            if (tx1 - tx0 + 1) * (ty1 - ty0 + 1) <= MAX_TILES or z == 0:
                break
            z -= 1
        self.zoom = z
        rows = []
        for ty in range(ty0, ty1 + 1):
            row = [self._tile(z, tx % 2 ** z, ty) for tx in range(tx0, tx1 + 1)]
            rows.append(np.concatenate(row, axis=1))
        self.mosaic, self.despiked = despike(np.concatenate(rows, axis=0),
                                             max(200.0, 2.0 * self.resolution(latc)))
        self.tx0, self.ty0 = tx0, ty0

    def _tile(self, z, x, y):
        import numpy as np
        path = os.path.join(self.cache, str(z), str(x), "%d.png" % y)
        if os.path.exists(path):
            with open(path, "rb") as f:
                data = f.read()
            try:
                rgb = decode_png(data)
            except TerrainError as e:
                raise TerrainError("cached tile %s is damaged (%s); delete it and run again"
                                   % (path, e))
        else:
            data = self.fetch(TILE_URL.format(z=z, x=x, y=y))
            rgb = decode_png(data)          # a damaged download is not cached
            os.makedirs(os.path.dirname(path), exist_ok=True)
            tmp = "%s.%d.tmp" % (path, os.getpid())
            with open(tmp, "wb") as f:
                f.write(data)
            os.replace(tmp, path)
        rgb = rgb.astype(np.float32)
        return rgb[..., 0] * 256.0 + rgb[..., 1] + rgb[..., 2] / 256.0 - 32768.0

    def sample(self, lat, lon):
        fx, fy = _tile_xy(lat, self._unwrap(lon), self.zoom)
        px = (fx - self.tx0) * 256.0 - 0.5
        py = (fy - self.ty0) * 256.0 - 0.5
        return _bilinear(self.mosaic, px, py)


class AnalyticDEM:
    """A function of (lat, lon) -> elevation, for tests."""

    def __init__(self, func, res=5.0, name="analytic"):
        self.func, self.res, self.name = func, res, name

    def resolution(self, lat):
        return self.res

    def sample(self, lat, lon):
        return self.func(lat, lon)


def despike(h, threshold):
    """Replace pixels more than threshold [m] from the median of their 3 x 3 neighbourhood by that
    median (the tiles carry isolated spikes and holes).  Returns (field, pixels replaced)."""
    import numpy as np
    pad = np.pad(h, 1, mode="edge")
    ny, nx = h.shape
    stack = np.stack([pad[j:j + ny, i:i + nx] for j in range(3) for i in range(3)])
    med = np.median(stack, axis=0)
    bad = np.abs(h - med) > threshold
    out = np.where(bad, med, h)
    return out, int(bad.sum())


def _http_get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "erf-domain-advisor"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.read()
    except (OSError, http.client.HTTPException) as e:
        raise TerrainError("could not download %s (%s); use --dem none or an ESRI ASCII grid"
                           % (url, e))


def _tile_xy(lat, lon, z):
    import numpy as np
    lat = np.clip(np.asarray(lat, dtype=float), -85.0511, 85.0511)
    n = 2.0 ** z
    x = (np.asarray(lon, dtype=float) + 180.0) / 360.0 * n
    y = (1.0 - np.log(np.tan(np.radians(lat)) + 1.0 / np.cos(np.radians(lat))) / math.pi) / 2.0 * n
    return x, y


def _bilinear(z, fx, fy):
    """Bilinear sample of z[row=y, col=x] at fractional indices; clamps at the edges."""
    import numpy as np
    ny, nx = z.shape
    fx = np.clip(fx, 0, nx - 1.000001)
    fy = np.clip(fy, 0, ny - 1.000001)
    i, j = np.floor(fx).astype(int), np.floor(fy).astype(int)
    tx, ty = fx - i, fy - j
    return ((1 - tx) * (1 - ty) * z[j, i] + tx * (1 - ty) * z[j, i + 1] +
            (1 - tx) * ty * z[j + 1, i] + tx * ty * z[j + 1, i + 1])


# --------------------------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------------------------

def level_terrain(frame, lv, dem, latc):
    """Cell-mean terrain on level lv's grid from k x k samples per cell.  Returns (h, water
    fraction, spacing): a level of more than MAX_SAMPLES cells is sampled every st cells (one
    sample per sampled cell, spacing st dx), which bounds the memory."""
    import numpy as np
    if hasattr(dem, "prepare"):
        # choose the source's zoom for this level first; k below depends on it
        cx = lv.x0 + (np.array([0.0, 0.5, 1.0])[None, :]) * (lv.x1 - lv.x0)
        cy = lv.y0 + (np.array([0.0, 0.5, 1.0])[:, None]) * (lv.y1 - lv.y0)
        clat, clon = frame.to_latlon_arrays(*np.broadcast_arrays(cx, cy))
        dem.prepare(clat, clon, lv.dx)
    st = max(1, int(math.ceil(math.sqrt(lv.nx * lv.ny / MAX_SAMPLES))))
    nx, ny, dx = lv.nx // st, lv.ny // st, lv.dx * st
    k = 1 if st > 1 else int(min(4, max(1, math.ceil(dx / max(dem.resolution(latc), 1e-3)))))
    while nx * ny * k * k > MAX_SAMPLES and k > 1:
        k -= 1
    sub = (np.arange(k) + 0.5) / k
    xs = lv.x0 + (np.arange(nx)[:, None] + sub[None, :]).ravel() * dx
    ys = lv.y0 + (np.arange(ny)[:, None] + sub[None, :]).ravel() * dx
    X, Y = np.meshgrid(xs, ys)
    lat, lon = frame.to_latlon_arrays(X, Y)
    if hasattr(dem, "covers") and not dem.covers(lat, lon):
        raise TerrainError("the elevation file %s does not cover level %d" % (dem.name, lv.lev))
    h = dem.sample(lat, lon).reshape(ny, k, nx, k)
    water = (h <= 0.0).mean(axis=(1, 3))
    # WPS terrain over the sea is 0 m: bathymetry would add slopes no grid sees.
    return np.maximum(h, 0.0).mean(axis=(1, 3)), water, dx


def cell_slope(h, dx):
    """Slope magnitude per cell from one-cell drops, as the grid sees them: in each direction the
    larger of the drops to the two neighbours over dx, then the magnitude of the two."""
    import numpy as np
    gx = np.zeros_like(h, dtype=float)
    gy = np.zeros_like(h, dtype=float)
    if h.shape[1] > 1:
        d = np.abs(np.diff(h, axis=1)) / dx
        gx[:, :-1] = d
        gx[:, 1:] = np.maximum(gx[:, 1:], d)
    if h.shape[0] > 1:
        d = np.abs(np.diff(h, axis=0)) / dx
        gy[:-1, :] = d
        gy[1:, :] = np.maximum(gy[1:, :], d)
    return np.hypot(gx, gy)


def slope_stats(h, dx, caution, limit, water=None):
    import numpy as np
    s = cell_slope(h, dx)
    finite = s[np.isfinite(s)]
    if finite.size == 0:
        raise TerrainError("no valid elevation values")
    return {
        "max_slope": float(finite.max()),
        "p99_slope": float(np.percentile(finite, 99)),
        "frac_caution": float((finite > caution).mean()),
        "frac_limit": float((finite > limit).mean()),
        "relief": float(np.nanmax(h) - np.nanmin(h)),
        "h_min": float(np.nanmin(h)),
        "h_max": float(np.nanmax(h)),
        "water_frac": float(np.mean(water if water is not None else (h <= 0.0))),
    }


def box_filter(h, w):
    """Mean over a w x w window (w odd); values beyond the edges repeat the edge value."""
    import numpy as np
    if w <= 1:
        return h.copy()
    r = w // 2
    pad = np.pad(h, r, mode="edge")
    c = np.cumsum(np.cumsum(pad, axis=0), axis=1)
    c = np.pad(c, ((1, 0), (1, 0)))
    ny, nx = h.shape
    return (c[w:w + ny, w:w + nx] - c[:ny, w:w + nx] - c[w:w + ny, :nx] + c[:ny, :nx]) / (w * w)


def boundary_mismatch(h0, dx0, width):
    """Largest |h0 - h0 averaged over 3 km| over the outer `width` cells of level 0: how far the
    level-0 terrain is from HRRR-like terrain where the boundary data enter."""
    import numpy as np
    w = max(1, int(round(3000.0 / dx0)))
    w += 1 - w % 2
    d = np.abs(h0 - box_filter(h0, w))
    mask = np.zeros_like(d, dtype=bool)
    mask[:width, :] = mask[-width:, :] = mask[:, :width] = mask[:, -width:] = True
    return float(np.nanmax(d[mask]))


def check(plan, dem, caution=0.3, limit=0.6, mismatch=150.0, keep=None):
    """Per-level statistics and warnings.  Returns (stats per level, warnings).  With a dict
    keep, keep[lev] = (terrain, water fraction, spacing) for --plot."""
    warns, stats = [], []
    latc = plan.settings.lat
    for lv in plan.levels:
        h, water, dxs = level_terrain(plan.frame, lv, dem, latc)
        if keep is not None:
            keep[lv.lev] = (h, water, dxs)
        st = slope_stats(h, dxs, caution, limit, water)
        if dxs > lv.dx:
            warns.append("Level %d has %d x %d cells: its terrain was sampled every %d cells, so "
                         "its slopes are underestimated." % (lv.lev, lv.nx, lv.ny, round(dxs / lv.dx)))
        # over terrain of height h the basic terrain-following grid squeezes the first cell to
        # dz0 (top - h) / top, so the slope factor is largest over the highest ground
        depth = plan.top - st["h_max"]
        room = plan.zi + damping_depth(plan.top)
        if depth < room and plan.top >= room:   # a top below room is warned about in the plan
            warns.append("Level %d: the highest terrain (%.0f m) leaves %.0f m to the model top "
                         "(%.0f m), less than the boundary layer (%.0f m) plus the damping layer "
                         "(%.0f m); raise --top." % (lv.lev, st["h_max"], depth, plan.top,
                                                    plan.zi, damping_depth(plan.top)))
        dz1 = plan.dz0 * max(depth, 0.0) / plan.top
        st["alpha"] = st["max_slope"] * lv.dx / dz1 if dz1 > 0 else float("inf")
        if st["alpha"] > 1.0:
            warns.append("Level %d: the terrain slope factor alpha = slope dx / dz1 reaches %.1f "
                         "(slope %.2f, dx %g m, first cell dz1 = %.3g m over the highest ground). "
                         "Above 1 the explicit terrain-metric diffusion limits dt to about "
                         "dz^2 / (4 K slope^2); ERF reports alpha in its log when "
                         "erf.diffusive_dt_check = true. erf.implicit_terrain_metric = true (with "
                         "erf.vert_implicit_fac > 0) takes that term into the implicit vertical "
                         "solve." % (lv.lev, st["alpha"], st["max_slope"], lv.dx, dz1))
        st["despiked"] = getattr(dem, "despiked", 0)
        st["lev"], st["dx"], st["dem_res"] = lv.lev, lv.dx, dem.resolution(latc)
        if lv.lev == 0:
            rim = max(1, int(math.ceil(plan.settings.real_width * lv.dx / dxs)))
            st["boundary_mismatch"] = boundary_mismatch(h, dxs, rim)
            if st["boundary_mismatch"] > mismatch:
                warns.append("Along the outer relaxation zone the level-0 terrain differs from "
                             "3 km-averaged (HRRR-like) terrain by up to %.0f m: the boundary "
                             "columns will be extrapolated. Move the domain edges off the steep "
                             "terrain, or blend the terrain toward the coarse terrain near the "
                             "edges in WPS." % st["boundary_mismatch"])
            if 0.02 < st["water_frac"] < 0.98:
                warns.append("The outer domain is %.0f%% water (elevation at or below 0 m; inland "
                             "basins below sea level count too): check "
                             "the land mask and the sea-surface temperature in the met_em files."
                             % (100 * st["water_frac"]))
        if st["max_slope"] > limit:
            warns.append("Level %d (dx %g m): the steepest slope is %.2f (%.0f deg), %.1f%% of "
                         "the cells above %.2f. A terrain-following grid may fail there: smooth the "
                         "terrain in WPS (smooth_option, smooth_passes), keep this level coarser, "
                         "or move its box. Metgrid runs need erf.terrain_type = StaticFittedMesh, so "
                         "immersed terrain is not an option here."
                         % (lv.lev, lv.dx, st["max_slope"], math.degrees(math.atan(st["max_slope"])),
                            100 * st["frac_limit"], limit))
        elif st["max_slope"] > caution:
            warns.append("Level %d (dx %g m): the steepest slope is %.2f (%.0f deg); above %.2f, "
                         "check the first steps for noise near the steep cells."
                         % (lv.lev, lv.dx, st["max_slope"], math.degrees(math.atan(st["max_slope"])),
                            caution))
        if lv.dx < 1000.0 and lv.lev > 0 and lv.lev == len(plan.levels) - 1:
            if st["dem_res"] > lv.dx:
                warns.append("The elevation data (%.0f m) are coarser than the finest level (%g m); "
                             "its slopes are underestimated." % (st["dem_res"], lv.dx))
        stats.append(st)
    return stats, warns


def make_dem(spec):
    if spec in (None, "none"):
        return None
    if spec == "terrarium":
        return TerrariumDEM()
    if not os.path.exists(spec):
        raise TerrainError("--dem %s: no such file" % spec)
    return AsciiGridDEM(spec)
