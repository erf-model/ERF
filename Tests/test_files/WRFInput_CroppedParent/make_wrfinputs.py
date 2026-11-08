#!/usr/bin/env python3
"""Write small, self-contained WRF inputs with a cropped middle ERF level.

d02 covers level-1 cells 18..53.  The test inputs crop ERF level 1 to 24..47.
d03 begins at level-2 cell 78 = 3 * (18 + 8); using the cropped ERF patch
origin instead would incorrectly place it at 96 = 3 * (24 + 8).
"""

from pathlib import Path

import netCDF4 as nc
import numpy as np


DATE = "2026-01-01_00:00:00"
NZ = 8
GRAV = 9.81


def write_file(path, nx, ny, ratio, parent_start, global_origin, total_ratio, dx):
    path = Path(path)
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
            ("I_PARENT_START", parent_start + 1),
            ("J_PARENT_START", parent_start + 1),
            ("PARENT_GRID_RATIO", ratio),
            ("USE_THETA_M", 0),
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

        gi = global_origin + np.arange(nx, dtype=np.float32)[None, :]
        gj = global_origin + np.arange(ny, dtype=np.float32)[:, None]
        # The level-2 theta pattern is checked in the plotfile at global indices.
        theta_pert = 0.01 * gi + 0.001 * gj
        put("ALT", xyz, 1.0)
        put("THM", xyz, theta_pert[None, None, :, :])
        put("U", ("Time", "bottom_top", "south_north", "west_east_stag"), 0.0)
        put("V", ("Time", "bottom_top", "south_north_stag", "west_east"), 0.0)
        put("W", zface, 0.0)
        put("PH", zface, 0.0)
        height = (GRAV * 100.0 * np.arange(NZ + 1, dtype=np.float32))[None, :, None, None]
        put("PHB", zface, height)
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
        put("XLAT", xy, 30.0 + 0.001 * (gj + 0.5) / total_ratio)
        put("XLONG", xy, -100.0 + 0.001 * (gi + 0.5) / total_ratio)
        put("IVGTYP", xy, 16, "i4")
        put("ISLTYP", xy, 14, "i4")


def write_low_file(path, nx, ny):
    with nc.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("Time", 1)
        ds.createDimension("DateStrLen", len(DATE))
        ds.createDimension("south_north", ny)
        ds.createDimension("west_east", nx)
        ds.createVariable("Times", "S1", ("Time", "DateStrLen"))[:] = \
            np.array(list(DATE), dtype="S1")
        ds.createVariable("SST", "f4", ("Time", "south_north", "west_east"))[:] = 300.0


def main():
    write_file("wrfinput_d01", 48, 48, 1, 0, 0, 1, 1000.0)
    write_file("wrfinput_d02", 36, 36, 3, 6, 18, 3, 1000.0 / 3.0)
    write_file("wrfinput_d03", 18, 18, 3, 8, 78, 9, 1000.0 / 9.0)
    write_low_file("wrflowinp_d01", 48, 48)
    write_low_file("wrflowinp_d02", 36, 36)
    write_low_file("wrflowinp_d03", 18, 18)


if __name__ == "__main__":
    main()
