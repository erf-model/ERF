"""Lambert conformal projection on a sphere, the HRRR CONUS grid, and the WPS domain frame.

The HRRR CONUS grid (GRIB2 grid definition): Lambert conformal, both standard parallels at
38.5 N, central meridian 97.5 W, a sphere of radius 6371229 m, 1799 x 1059 points 3 km apart,
first point at 21.138123 N, 122.719528 W.
With a radius of 6370000 m the same first point puts the grid centre exactly on 38.5 N,
97.5 W; with 6371229 m it is about 0.6 km (about 0.2 HRRR cells) away.  The HRRR edge margin is therefore known
to within about one HRRR cell, well inside the default --hrrr-buffer of 5 cells.

WPS geogrid uses a sphere of radius 6370000 m.  The advisor puts the WPS domain on the same
standard parallels and central meridian as HRRR, so the two grids are not rotated against each
other; outside the HRRR domain it uses the point's own latitude and longitude instead.
"""

import math

HRRR_RADIUS = 6371229.0
WPS_RADIUS = 6370000.0
HRRR_TRUELAT = 38.5
HRRR_STAND_LON = -97.5
HRRR_NX = 1799
HRRR_NY = 1059
HRRR_DX = 3000.0
HRRR_LAT1 = 21.138123
HRRR_LON1 = -122.719528


class Lambert:
    """Tangent Lambert conformal projection: x east, y north, metres, origin at (lat0, lon0)."""

    def __init__(self, truelat, stand_lon, lat0=None, radius=WPS_RADIUS):
        if not (-89.0 <= truelat <= 89.0) or abs(truelat) < 1.0:
            raise ValueError("Lambert standard parallel must lie between 1 and 89 degrees "
                             "from the equator, got %g" % truelat)
        self.radius = radius
        self.lon0 = stand_lon
        p1 = math.radians(truelat)
        self.n = math.sin(p1)
        self.F = math.cos(p1) * math.tan(math.pi / 4 + p1 / 2) ** self.n / self.n
        self.rho0 = self._rho(truelat if lat0 is None else lat0)

    def _rho(self, lat):
        phi = math.radians(lat)
        return self.radius * self.F / math.tan(math.pi / 4 + phi / 2) ** self.n

    def forward(self, lat, lon):
        dlon = (lon - self.lon0 + 180.0) % 360.0 - 180.0
        theta = self.n * math.radians(dlon)
        rho = self._rho(lat)
        return rho * math.sin(theta), self.rho0 - rho * math.cos(theta)

    def inverse(self, x, y):
        dy = self.rho0 - y
        rho = math.copysign(math.hypot(x, dy), self.n)
        theta = math.atan2(math.copysign(1.0, self.n) * x, math.copysign(1.0, self.n) * dy)
        lat = math.degrees(2 * math.atan((self.radius * self.F / rho) ** (1.0 / self.n)) - math.pi / 2)
        lon = self.lon0 + math.degrees(theta / self.n)
        return lat, (lon + 180.0) % 360.0 - 180.0

    def inverse_arrays(self, x, y):
        """numpy version of inverse()."""
        import numpy as np
        dy = self.rho0 - y
        s = 1.0 if self.n > 0 else -1.0
        rho = s * np.hypot(x, dy)
        theta = np.arctan2(s * x, s * dy)
        lat = np.degrees(2 * np.arctan((self.radius * self.F / rho) ** (1.0 / self.n)) - np.pi / 2)
        lon = self.lon0 + np.degrees(theta / self.n)
        return lat, (lon + 180.0) % 360.0 - 180.0

    def map_factor(self, lat):
        phi = math.radians(lat)
        return self._rho(lat) * self.n / (self.radius * math.cos(phi))


HRRR = Lambert(HRRR_TRUELAT, HRRR_STAND_LON, radius=HRRR_RADIUS)
_HRRR_X1, _HRRR_Y1 = HRRR.forward(HRRR_LAT1, HRRR_LON1)


def hrrr_index(lat, lon):
    """Fractional (i, j) of a point on the HRRR grid, 0-based from the first point."""
    x, y = HRRR.forward(lat, lon)
    return (x - _HRRR_X1) / HRRR_DX, (y - _HRRR_Y1) / HRRR_DX


def hrrr_margin_cells(lat, lon):
    """Distance from the point to the nearest HRRR grid edge in HRRR cells (negative outside)."""
    i, j = hrrr_index(lat, lon)
    return min(i, j, HRRR_NX - 1 - i, HRRR_NY - 1 - j)


class DomainFrame:
    """The WPS d01 / ERF level-0 frame: ERF coordinates (x, y) in metres from the lower-left
    corner of a domain of size lx by ly centred on (ref_lat, ref_lon)."""

    def __init__(self, ref_lat, ref_lon, lx, ly, truelat, stand_lon):
        self.ref_lat, self.ref_lon = ref_lat, ref_lon
        self.lx, self.ly = lx, ly
        self.truelat, self.stand_lon = truelat, stand_lon
        self.proj = Lambert(truelat, stand_lon, radius=WPS_RADIUS)
        self.xc, self.yc = self.proj.forward(ref_lat, ref_lon)

    def to_erf(self, lat, lon):
        x, y = self.proj.forward(lat, lon)
        return x - self.xc + self.lx / 2, y - self.yc + self.ly / 2

    def to_latlon(self, x, y):
        return self.proj.inverse(x - self.lx / 2 + self.xc, y - self.ly / 2 + self.yc)

    def to_latlon_arrays(self, x, y):
        return self.proj.inverse_arrays(x - self.lx / 2 + self.xc, y - self.ly / 2 + self.yc)


def choose_projection(lat, lon):
    """HRRR's standard parallel and meridian inside the HRRR domain, else the point's own."""
    if hrrr_margin_cells(lat, lon) >= 0:
        return HRRR_TRUELAT, HRRR_STAND_LON, True
    truelat = lat if abs(lat) >= 1.0 else math.copysign(1.0, lat if lat != 0 else 1.0)
    return truelat, lon, False
