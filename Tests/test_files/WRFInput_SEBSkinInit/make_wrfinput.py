#!/usr/bin/env python3
"""Write a small wrfinput file whose skin temperature TSK and deep soil temperature TMN
vary from column to column.

The state is at rest and uniform. TSK = 280 + 0.5 i + 0.25 j [K] and
TMN = 275 + 0.25 i + 0.5 j [K] (i, j column indices from zero), and LANDMASK is 0 (water)
where i >= 6 and j <= 2, a block that lies in both boxes of the deck. check_skin.py
imports the same formulas.

  --all-land           no water column
  --land VALUE         TSK of the land column (8, 5) is VALUE
  --land-west VALUE    TSK of the land column (2, 5), in the other box, is VALUE
  --water VALUE        TSK of the water column (8, 1) is VALUE
  --nest               also write wrfinput_d02 and wrflowinp_d02: 8 x 8 columns refined by 2 over the level 0
                       columns 2..5 x 2..5, with TSK = 290 + 0.125 I + 0.0625 J [K] and
                       TMN = 284 + 0.0625 I + 0.125 J [K] (I, J the fine column indices of
                       the whole domain, 4..11) and water in the fine columns (10, 5) and
                       (11, 5)
  --no-tmn             leave TMN out of the files
  --tmn-land VALUE     TMN of the land column (8, 5) is VALUE
  --tmn-land-west VALUE  TMN of the land column (2, 5), in the other box, is VALUE
  --tmn-water VALUE    TMN of the water column (8, 1) is VALUE
"""

import argparse

import netCDF4 as nc
import numpy as np

DATE = "2026-01-01_00:00:00"
NX, NY, NZ = 12, 8, 8
GRAV = 9.81
RAD_T_SFC = 300.0        # erf.rad_t_sfc of the deck [K]
T_DEEP_DEFAULT = 290.0   # erf.radiation.seb_t_deep_default of the deck [K]
SST = 285.0              # the sea surface temperature of the file [K]
LAND_COLUMN = (8, 5)     # (i, j): the first column of the deck's second box
WEST_LAND_COLUMN = (2, 5)  # a land column of the deck's first box
WATER_COLUMN = (8, 1)


def tsk(i, j):
    """The skin temperature of the file [K] at column (i, j)."""
    return 280.0 + 0.5 * i + 0.25 * j


def tmn(i, j):
    """The deep soil temperature of the file [K] at column (i, j)."""
    return 275.0 + 0.25 * i + 0.5 * j


def is_water(i, j, all_land=False):
    return (not all_land) and i >= 6 and j <= 2


# The nest of --nest: level 0 columns NEST_LO..NEST_LO+3 in x and y, refined by NEST_RATIO
NEST_LO, NEST_RATIO, NEST_N = 2, 2, 8


def nest_tsk(i, j):
    """The skin temperature of wrfinput_d02 [K] at the fine column (i, j) of the whole domain."""
    return 290.0 + 0.125 * i + 0.0625 * j


def nest_tmn(i, j):
    """The deep soil temperature of wrfinput_d02 [K] at the fine column (i, j) of the whole domain."""
    return 284.0 + 0.0625 * i + 0.125 * j


def nest_is_water(i, j):
    return j == 5 and i >= 10


def write_file(path, nx, ny, ratio, parent_start, dx, skin, landmask, deep=None):
    """skin, deep: (ny, nx) arrays; None leaves TSK or TMN out."""
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
            ("I_PARENT_START", parent_start + 1), ("J_PARENT_START", parent_start + 1),
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
        put("PHB", zface, (GRAV * 100.0 * np.arange(NZ + 1, dtype=np.float32))[None, :, None, None])
        put("PB", xyz, 100000.0)
        put("P", xyz, 0.0)
        for name, value in (("PSFC", 100000.0), ("MUB", 10000.0), ("MAPFAC_M", 1.0), ("SST", SST)):
            put(name, xy, value)
        if skin is not None:
            put("TSK", xy, skin[None, ...])
        if deep is not None:
            put("TMN", xy, deep[None, ...])
        put("LANDMASK", xy, landmask[None, ...])
        put("MAPFAC_U", ("Time", "south_north", "west_east_stag"), 1.0)
        put("MAPFAC_V", ("Time", "south_north_stag", "west_east"), 1.0)
        put("C1H", ("Time", "bottom_top"), 1.0)
        put("C2H", ("Time", "bottom_top"), 0.0)
        put("RDNW", ("Time", "bottom_top"), 1.0 / NZ)
        put("XLAT", xy, 30.0)
        put("XLONG", xy, -100.0)
        put("IVGTYP", xy, 16, "i4")
        put("ISLTYP", xy, 14, "i4")


def write_low_file(path, nx, ny):
    """The lower boundary file: one time, a uniform sea surface temperature."""
    with nc.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("Time", 1)
        ds.createDimension("DateStrLen", len(DATE))
        ds.createDimension("south_north", ny)
        ds.createDimension("west_east", nx)
        ds.createVariable("Times", "S1", ("Time", "DateStrLen"))[:] = np.array(list(DATE), dtype="S1")
        ds.createVariable("SST", "f4", ("Time", "south_north", "west_east"))[:] = SST


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-land", action="store_true")
    parser.add_argument("--land", type=float, default=None)
    parser.add_argument("--land-west", type=float, default=None)
    parser.add_argument("--water", type=float, default=None)
    parser.add_argument("--nest", action="store_true")
    parser.add_argument("--no-tmn", action="store_true")
    parser.add_argument("--tmn-land", type=float, default=None)
    parser.add_argument("--tmn-land-west", type=float, default=None)
    parser.add_argument("--tmn-water", type=float, default=None)
    args = parser.parse_args()

    i = np.arange(NX, dtype=np.float64)[None, :] + np.zeros((NY, 1))
    j = np.arange(NY, dtype=np.float64)[:, None] + np.zeros((1, NX))
    skin = tsk(i, j)
    landmask = np.ones((NY, NX))
    if not args.all_land:
        landmask[(i >= 6) & (j <= 2)] = 0.0
    if args.land is not None:
        skin[LAND_COLUMN[1], LAND_COLUMN[0]] = args.land
    if args.land_west is not None:
        skin[WEST_LAND_COLUMN[1], WEST_LAND_COLUMN[0]] = args.land_west
    if args.water is not None:
        skin[WATER_COLUMN[1], WATER_COLUMN[0]] = args.water

    deep = None if args.no_tmn else tmn(i, j)
    if deep is not None:
        for column, value in ((LAND_COLUMN, args.tmn_land), (WEST_LAND_COLUMN, args.tmn_land_west),
                              (WATER_COLUMN, args.tmn_water)):
            if value is not None:
                deep[column[1], column[0]] = value

    write_file("wrfinput_d01", NX, NY, 1, 0, 1000.0, skin, landmask, deep)

    if args.nest:
        fi = NEST_RATIO * NEST_LO + np.arange(NEST_N, dtype=np.float64)[None, :] + np.zeros((NEST_N, 1))
        fj = NEST_RATIO * NEST_LO + np.arange(NEST_N, dtype=np.float64)[:, None] + np.zeros((1, NEST_N))
        fine_mask = np.ones((NEST_N, NEST_N))
        fine_mask[(fj == 5) & (fi >= 10)] = 0.0
        write_file("wrfinput_d02", NEST_N, NEST_N, NEST_RATIO, NEST_LO, 1000.0 / NEST_RATIO,
                   nest_tsk(fi, fj), fine_mask, None if args.no_tmn else nest_tmn(fi, fj))

    write_low_file("wrflowinp_d01", NX, NY)
    if args.nest:
        write_low_file("wrflowinp_d02", NEST_N, NEST_N)


if __name__ == "__main__":
    main()
