#!/usr/bin/env python3
"""Write a small, self-contained wrfinput whose PSFC, SST and TSK vary in space.

PSFC = 100000 + 150 i - 100 j [Pa] (exact in single precision), SST = 295 + 0.2 i +
0.1 j [K] and TSK = SST + 1 + 0.05 i, so that converting a cell with the pressure of
another cell, or of another file, shows.  wrfinput_d01b is the same with PSFC 5000 Pa
higher: listed second in erf.nc_init_file_0, it must not be the file the level-0
surface pressure is recovered from.
"""

from pathlib import Path

import netCDF4 as nc
import numpy as np


DATE = "2026-01-01_00:00:00"
NX = 16
NY = 16
NZ = 8
DX = 1000.0
GRAV = 9.81


def write_file(path, psfc_offset):
    with nc.Dataset(Path(path), "w", format="NETCDF4_CLASSIC") as ds:
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
            ("I_PARENT_START", 1),
            ("J_PARENT_START", 1),
            ("PARENT_GRID_RATIO", 1),
            ("USE_THETA_M", 0),
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
        sst = 295.0 + 0.2 * i + 0.1 * j
        put("PSFC", xy, (100000.0 + psfc_offset + 150.0 * i - 100.0 * j)[None])
        put("SST", xy, sst[None])
        put("TSK", xy, (sst + 1.0 + 0.05 * i)[None])

        put("ALT", xyz, 1.0)
        put("THM", xyz, 0.0)
        put("QVAPOR", xyz, 0.0)
        put("U", ("Time", "bottom_top", "south_north", "west_east_stag"), 0.0)
        put("V", ("Time", "bottom_top", "south_north_stag", "west_east"), 0.0)
        put("W", zface, 0.0)
        put("PH", zface, 0.0)
        height = (GRAV * 100.0 * np.arange(NZ + 1, dtype=np.float32))[None, :, None, None]
        put("PHB", zface, height)
        put("PB", xyz, 100000.0)
        put("P", xyz, 0.0)
        for name, value in (("MUB", 10000.0), ("MAPFAC_M", 1.0), ("LANDMASK", 0.0)):
            put(name, xy, value)
        put("MAPFAC_U", ("Time", "south_north", "west_east_stag"), 1.0)
        put("MAPFAC_V", ("Time", "south_north_stag", "west_east"), 1.0)
        put("C1H", ("Time", "bottom_top"), 1.0)
        put("C2H", ("Time", "bottom_top"), 0.0)
        put("RDNW", ("Time", "bottom_top"), 1.0 / NZ)
        put("XLAT", xy, 30.0 + 0.001 * (j + 0.5))
        put("XLONG", xy, -100.0 + 0.001 * (i + 0.5))
        put("IVGTYP", xy, 16, "i4")
        put("ISLTYP", xy, 14, "i4")


def main():
    write_file("wrfinput_d01", 0.0)
    write_file("wrfinput_d01b", 5000.0)


if __name__ == "__main__":
    main()
