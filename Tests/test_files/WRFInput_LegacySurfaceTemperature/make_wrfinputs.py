#!/usr/bin/env python3
"""Write a small, self-contained wrfinput whose PSFC, SST and TSK vary in space.

PSFC = 100000 + 150 i - 100 j [Pa] (exact in single precision), SST = 295 + 0.2 i +
0.1 j [K] and TSK = SST + 1 + 0.05 i, so that converting a cell with the pressure of
another cell, or of another file, shows.  wrfinput_d01b is the same with PSFC 5000 Pa
higher: listed second in erf.nc_init_file_0, it must not be the file the level-0
surface pressure is recovered from.

wrfinput_mixed is wrfinput_d01 with land (LANDMASK 1) west of i = LAND_COLUMNS, where the
SST is the placeholder WRF writes over land, zero; wrfinput_mixed_badpsfc is that with
PSFC -1 Pa at the land cell BAD_PSFC_CELL (i, j).

nest_d01 (PSFC 2000 Pa lower) has two nests, read at level 1 by nest_files.i: nest_d02
(5000 Pa lower) over its cells 2..5, and nest_d03 (10000 Pa lower) over 10..13, which
starts 10 s later, so a run that stops before then never reads it.  i and j are the
indices of each file.
"""

from pathlib import Path

import netCDF4 as nc
import numpy as np


DATE = "2026-01-01_00:00:00"
NX = 16
NZ = 8
DX = 1000.0
GRAV = 9.81
LAND_COLUMNS = 8
BAD_PSFC_CELL = (3, 9)
NEST_RATIO = 3
NEST_CELLS = 12


def write_file(path, psfc_offset, land_columns=0, bad_psfc_cell=None, nx=NX, ratio=1,
               parent_start=0, date=DATE):
    ny = nx
    with nc.Dataset(Path(path), "w", format="NETCDF4_CLASSIC") as ds:
        for name, size in (
            ("Time", 1), ("DateStrLen", len(date)),
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
        ds.setncattr("DX", np.float32(DX / ratio))
        ds.setncattr("DY", np.float32(DX / ratio))
        ds.setncattr("START_DATE", date)
        ds.setncattr("SIMULATION_START_DATE", date)

        xy = ("Time", "south_north", "west_east")
        xyz = ("Time", "bottom_top", "south_north", "west_east")
        zface = ("Time", "bottom_top_stag", "south_north", "west_east")

        def put(name, dims, data, dtype="f4"):
            var = ds.createVariable(name, dtype, dims)
            var[:] = np.broadcast_to(data, var.shape)

        put("Times", ("Time", "DateStrLen"), np.array(list(date), dtype="S1"), "S1")

        i = np.arange(nx, dtype=np.float64)[None, :]
        j = np.arange(ny, dtype=np.float64)[:, None]
        sst = 295.0 + 0.2 * i + 0.1 * j
        land = np.broadcast_to(i < land_columns, (ny, nx))
        psfc = np.array(np.broadcast_to(100000.0 + psfc_offset + 150.0 * i - 100.0 * j,
                                        (ny, nx)))
        if bad_psfc_cell is not None:
            psfc[bad_psfc_cell[1], bad_psfc_cell[0]] = -1.0
        put("PSFC", xy, psfc[None])
        put("SST", xy, np.where(land, 0.0, sst)[None])
        put("TSK", xy, (sst + 1.0 + 0.05 * i)[None])
        put("LANDMASK", xy, land.astype(np.float64)[None])

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
        for name, value in (("MUB", 10000.0), ("MAPFAC_M", 1.0)):
            put(name, xy, value)
        put("MAPFAC_U", ("Time", "south_north", "west_east_stag"), 1.0)
        put("MAPFAC_V", ("Time", "south_north_stag", "west_east"), 1.0)
        put("C1H", ("Time", "bottom_top"), 1.0)
        put("C2H", ("Time", "bottom_top"), 0.0)
        put("RDNW", ("Time", "bottom_top"), 1.0 / NZ)
        # At the indices of the outermost domain
        origin = ratio * parent_start
        put("XLAT", xy, 30.0 + 0.001 * (origin + j + 0.5) / ratio)
        put("XLONG", xy, -100.0 + 0.001 * (origin + i + 0.5) / ratio)
        put("IVGTYP", xy, 16, "i4")
        put("ISLTYP", xy, 14, "i4")


def main():
    write_file("wrfinput_d01", 0.0)
    write_file("wrfinput_d01b", 5000.0)
    write_file("wrfinput_mixed", 0.0, LAND_COLUMNS)
    write_file("wrfinput_mixed_badpsfc", 0.0, LAND_COLUMNS, BAD_PSFC_CELL)
    write_file("nest_d01", -2000.0)
    write_file("nest_d02", -5000.0, nx=NEST_CELLS, ratio=NEST_RATIO, parent_start=2)
    write_file("nest_d03", -10000.0, nx=NEST_CELLS, ratio=NEST_RATIO, parent_start=10,
               date="2026-01-01_00:00:10")


if __name__ == "__main__":
    main()
