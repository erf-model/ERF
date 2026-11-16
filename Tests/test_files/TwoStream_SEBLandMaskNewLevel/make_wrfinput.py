#!/usr/bin/env python3
"""Write the wrfinput of TwoStream_SEBLandMaskNewLevel: land in the west, water in the east.

A 16 x 16 x 8 grid of 1 km cells, flat, at rest under a top at 800 m, with theta = 300 K
and qv = 0.005. LANDMASK is 1 for i <= 7 and 0 for i >= 8, except the land column (8, 6):
a level 1 over i = 6..11, j = 4..9 then holds land, a shore that is not a straight line,
and water three columns from the nearest land.

wrfinput_d02 is a file for that level 1 (12 x 12 cells of 500 m, parent start (7, 5) in
WRF's 1-based count): the same state, and the LANDMASK of the parent column under each
cell. One leg of the test gives it to level 1 as erf.nc_init_file_1.
"""

import netCDF4 as nc
import numpy as np

DATE = "2026-06-21_12:00:00"
NX = NY = 16
NZ = 8
DX = 1000.0
ZTOP = 800.0
GRAV = 9.81


def land(i, j):
    """LANDMASK of column (i, j)."""
    return 1 if (i <= 7 or (i == 8 and j == 6)) else 0


def write_input(path, nx=NX, ny=NY, dx=DX, ratio=1, parent_start=(0, 0)):
    """Write a wrfinput file of nx x ny cells of dx; cell (i, j) lies over the parent
    column (parent_start[0] + i // ratio, parent_start[1] + j // ratio)."""
    with nc.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        for name, size in (
            ("Time", 1), ("DateStrLen", len(DATE)),
            ("west_east", nx), ("west_east_stag", nx + 1),
            ("south_north", ny), ("south_north_stag", ny + 1),
            ("bottom_top", NZ), ("bottom_top_stag", NZ + 1),
        ):
            ds.createDimension(name, size)
        for name, value in (
            ("WEST-EAST_GRID_DIMENSION", nx + 1),
            ("SOUTH-NORTH_GRID_DIMENSION", ny + 1),
            ("BOTTOM-TOP_GRID_DIMENSION", NZ + 1),
            ("I_PARENT_START", parent_start[0] + 1), ("J_PARENT_START", parent_start[1] + 1),
            ("PARENT_GRID_RATIO", ratio), ("USE_THETA_M", 0),
        ):
            ds.setncattr(name, np.int32(value))
        ds.setncattr("DX", np.float32(dx))
        ds.setncattr("DY", np.float32(dx))
        ds.setncattr("START_DATE", DATE)
        ds.setncattr("SIMULATION_START_DATE", DATE)

        xy = ("Time", "south_north", "west_east")
        xyz = ("Time", "bottom_top", "south_north", "west_east")
        zface = ("Time", "bottom_top_stag", "south_north", "west_east")

        def put(name, dims, data, dtype="f4"):
            var = ds.createVariable(name, dtype, dims)
            var[:] = np.broadcast_to(data, var.shape)

        put("Times", ("Time", "DateStrLen"), np.array(list(DATE), dtype="S1"), "S1")
        put("ALT", xyz, 1.0)
        put("THM", xyz, 0.0)
        put("U", ("Time", "bottom_top", "south_north", "west_east_stag"), 0.0)
        put("V", ("Time", "bottom_top", "south_north_stag", "west_east"), 0.0)
        put("W", zface, 0.0)
        put("PH", zface, 0.0)
        height = ZTOP * np.arange(NZ + 1, dtype=np.float64)[:, None, None] / NZ
        put("PHB", zface, GRAV * height[None, :, :, :])
        put("QVAPOR", xyz, 0.005)
        put("PB", xyz, 100000.0)
        put("P", xyz, 0.0)
        mask = np.array([[land(parent_start[0] + i // ratio, parent_start[1] + j // ratio)
                          for i in range(nx)] for j in range(ny)], dtype=np.float64)
        for name, value in (
            ("PSFC", 100000.0), ("MUB", 10000.0),
            ("MAPFAC_M", 1.0), ("SST", 290.0), ("TSK", 290.0),
            ("LANDMASK", mask),
        ):
            put(name, xy, value)
        put("MAPFAC_U", ("Time", "south_north", "west_east_stag"), 1.0)
        put("MAPFAC_V", ("Time", "south_north_stag", "west_east"), 1.0)
        put("C1H", ("Time", "bottom_top"), 1.0)
        put("C2H", ("Time", "bottom_top"), 0.0)
        put("RDNW", ("Time", "bottom_top"), 1.0 / NZ)
        put("XLAT", xy, 40.0)
        put("XLONG", xy, -105.0)
        put("IVGTYP", xy, 16, "i4")
        put("ISLTYP", xy, 14, "i4")


def write_low_file(path, nx=NX, ny=NY):
    with nc.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("Time", 1)
        ds.createDimension("DateStrLen", len(DATE))
        ds.createDimension("south_north", ny)
        ds.createDimension("west_east", nx)
        ds.createVariable("Times", "S1", ("Time", "DateStrLen"))[:] = \
            np.array(list(DATE), dtype="S1")
        ds.createVariable("SST", "f4", ("Time", "south_north", "west_east"))[:] = 290.0


if __name__ == "__main__":
    write_input("wrfinput_d01")
    write_low_file("wrflowinp_d01")
    # Level 1 over coarse columns i = 6..11, j = 4..9.
    write_input("wrfinput_d02", nx=12, ny=12, dx=DX / 2, ratio=2, parent_start=(6, 4))
    write_low_file("wrflowinp_d02", nx=12, ny=12)
