#!/usr/bin/env python3
"""Write a small, self-contained wrfinput on a rotated latitude/longitude grid.

The 32 x 32 x 8 grid of 1 km cells is at rest under a flat top at 800 m. It is
flat, with 100 m levels, or with --hill its ground rises 20 m per cell in i
(ground() at a cell centre) and its 8 levels are spaced evenly from the ground
to the top. Its latitude and longitude are linear in (i, j) but turned from the
grid axes, so a row of cells is not a line of constant latitude.
theta = 300 + 0.1 i + 0.001 j names the column a sample comes from, and
qv = 0.001 (k + 1) its level. check_line_latlon.py uses latlon() to place the
lines it asks for.
"""

import sys

import netCDF4 as nc
import numpy as np


DATE = "2026-01-01_00:00:00"
NX = NY = 32
NZ = 8
DX = 1000.0
ZTOP = 800.0
GRAV = 9.81


def latlon(fi, fj):
    """Latitude and longitude at fractional cell coordinates (fi, fj); a cell
    centre is at (i + 0.5, j + 0.5)."""
    return (40.0 + 0.009 * fj + 0.002 * fi,
            -105.0 + 0.0117 * fi - 0.002 * fj)


def theta(i, j):
    return 300.0 + 0.1 * i + 0.001 * j


def ground(fi, hill):
    """Height of the ground at fractional cell coordinate fi in x."""
    return 20.0 * fi if hill else 0.0 * fi


def write_input(path, hill=False):
    with nc.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        for name, size in (
            ("Time", 1), ("DateStrLen", len(DATE)),
            ("west_east", NX), ("west_east_stag", NX + 1),
            ("south_north", NY), ("south_north_stag", NY + 1),
            ("bottom_top", NZ), ("bottom_top_stag", NZ + 1),
        ):
            ds.createDimension(name, size)
        for name, value in (
            ("WEST-EAST_GRID_DIMENSION", NX + 1),
            ("SOUTH-NORTH_GRID_DIMENSION", NY + 1),
            ("BOTTOM-TOP_GRID_DIMENSION", NZ + 1),
            ("I_PARENT_START", 1), ("J_PARENT_START", 1),
            ("PARENT_GRID_RATIO", 1), ("USE_THETA_M", 0),
        ):
            ds.setncattr(name, np.int32(value))
        ds.setncattr("DX", np.float32(DX))
        ds.setncattr("DY", np.float32(DX))
        ds.setncattr("START_DATE", DATE)
        ds.setncattr("SIMULATION_START_DATE", DATE)

        xy = ("Time", "south_north", "west_east")
        xyz = ("Time", "bottom_top", "south_north", "west_east")
        zface = ("Time", "bottom_top_stag", "south_north", "west_east")

        def put(name, dims, data, dtype="f4"):
            var = ds.createVariable(name, dtype, dims)
            var[:] = np.broadcast_to(data, var.shape)

        put("Times", ("Time", "DateStrLen"), np.array(list(DATE), dtype="S1"), "S1")

        i = np.arange(NX, dtype=np.float64)[None, :]
        j = np.arange(NY, dtype=np.float64)[:, None]
        put("ALT", xyz, 1.0)
        put("THM", xyz, (theta(i, j) - 300.0)[None, None, :, :])
        put("U", ("Time", "bottom_top", "south_north", "west_east_stag"), 0.0)
        put("V", ("Time", "bottom_top", "south_north_stag", "west_east"), 0.0)
        put("W", zface, 0.0)
        put("PH", zface, 0.0)
        h = ground(i + 0.5, hill)[None, :, :]
        height = h + (ZTOP - h) * np.arange(NZ + 1, dtype=np.float64)[:, None, None] / NZ
        put("PHB", zface, GRAV * height[None, :, :, :])
        put("QVAPOR", xyz, 0.001 * (np.arange(NZ, dtype=np.float64) + 1.0)[None, :, None, None])
        put("PB", xyz, 100000.0)
        put("P", xyz, 0.0)
        for name, value in (
            ("PSFC", 100000.0), ("MUB", 10000.0),
            ("MAPFAC_M", 1.0), ("SST", 300.0), ("TSK", 300.0),
            ("LANDMASK", 0.0),
        ):
            put(name, xy, value)
        put("MAPFAC_U", ("Time", "south_north", "west_east_stag"), 1.0)
        put("MAPFAC_V", ("Time", "south_north_stag", "west_east"), 1.0)
        put("C1H", ("Time", "bottom_top"), 1.0)
        put("C2H", ("Time", "bottom_top"), 0.0)
        put("RDNW", ("Time", "bottom_top"), 1.0 / NZ)
        lat, lon = latlon(i + 0.5, j + 0.5)
        put("XLAT", xy, lat)
        put("XLONG", xy, lon)
        put("IVGTYP", xy, 16, "i4")
        put("ISLTYP", xy, 14, "i4")


def write_low_file(path):
    with nc.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("Time", 1)
        ds.createDimension("DateStrLen", len(DATE))
        ds.createDimension("south_north", NY)
        ds.createDimension("west_east", NX)
        ds.createVariable("Times", "S1", ("Time", "DateStrLen"))[:] = \
            np.array(list(DATE), dtype="S1")
        ds.createVariable("SST", "f4", ("Time", "south_north", "west_east"))[:] = 300.0


if __name__ == "__main__":
    write_input("wrfinput_d01", hill="--hill" in sys.argv[1:])
    write_low_file("wrflowinp_d01")
