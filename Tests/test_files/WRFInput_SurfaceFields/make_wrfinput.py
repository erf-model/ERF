#!/usr/bin/env python3
"""Write a small wrfinput file that carries the land and soil fields.

The state is at rest and uniform. Every land and soil field varies in both horizontal
directions (and with the soil layer); check_log.py imports the same formulas.

  --omit NAME [NAME ...]   leave these variables out of the file
  --bad NAME VALUE         one column of the field holds VALUE (may be repeated)
  --fill NAME VALUE        every column of the field holds VALUE (may be repeated)
  --nz N                   atmospheric levels (default 8); 3 is fewer than the soil layers
  --zero-layers            ZS and DZS are declared and left at zero
  --surface-zs             ZS starts at zero (a land model with its first level at the surface)
"""

import argparse

import netCDF4 as nc
import numpy as np

DATE = "2026-01-01_00:00:00"
NX, NY, NSOIL = 12, 8, 4
GRAV = 9.81
# The four Noah layers [m]
DZS = np.array([0.1, 0.3, 0.6, 1.0], dtype=np.float32)
ZS = np.array([0.05, 0.25, 0.7, 1.5], dtype=np.float32)
# Depths that start at the surface [m]
ZS_SURFACE = np.array([0.0, 0.1, 0.4, 1.0], dtype=np.float32)


def surface_fields():
    """name -> array of shape (NY, NX) or (NSOIL, NY, NX)."""
    i = np.arange(NX, dtype=np.float64)[None, :] + np.zeros((NY, 1))
    j = np.arange(NY, dtype=np.float64)[:, None] + np.zeros((1, NX))
    k = np.arange(NSOIL, dtype=np.float64)[:, None, None]
    seaice = np.zeros((NY, NX))
    seaice[5, 9] = 1.0
    return {
        "LAI": 0.5 + 0.25 * i + 0.03125 * j,
        "VEGFRA": 10.0 + 5.0 * j + 0.25 * i,
        "TMN": 280.0 + 0.5 * i + 0.0625 * j,
        "ZNT": 0.01 + 0.01 * j + 0.001 * i,
        "ALBEDO": 0.1 + 0.02 * i + 0.001 * j,
        "ALBBCK": 0.15 + 0.01 * i + 0.002 * j,
        "SEAICE": seaice,
        "SNOW": 2.0 * j + 0.125 * i,
        "SNOWH": 0.01 * j + 0.0005 * i,
        "TSLB": 290.0 - 2.0 * k + 0.5 * i[None, :, :] + 0.0625 * j[None, :, :],
        "SMOIS": 0.1 + 0.05 * k + 0.003 * j[None, :, :] + 0.0002 * i[None, :, :],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omit", nargs="*", default=[])
    parser.add_argument("--bad", nargs=2, action="append", default=[], metavar=("NAME", "VALUE"))
    parser.add_argument("--fill", nargs=2, action="append", default=[], metavar=("NAME", "VALUE"))
    parser.add_argument("--surface-zs", action="store_true")
    parser.add_argument("--nz", type=int, default=8)
    parser.add_argument("--zero-layers", action="store_true")
    args = parser.parse_args()
    NZ = args.nz

    with nc.Dataset("wrfinput_d01", "w", format="NETCDF4_CLASSIC") as ds:
        for name, size in (
            ("Time", 1), ("DateStrLen", len(DATE)),
            ("west_east", NX), ("west_east_stag", NX + 1),
            ("south_north", NY), ("south_north_stag", NY + 1),
            ("bottom_top", NZ), ("bottom_top_stag", NZ + 1),
            ("soil_layers_stag", NSOIL),
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
        ds.setncattr("DX", np.float32(1000.0))
        ds.setncattr("DY", np.float32(1000.0))
        ds.setncattr("START_DATE", DATE)
        ds.setncattr("SIMULATION_START_DATE", DATE)

        xy = ("Time", "south_north", "west_east")
        xyz = ("Time", "bottom_top", "south_north", "west_east")
        zface = ("Time", "bottom_top_stag", "south_north", "west_east")
        soil = ("Time", "soil_layers_stag", "south_north", "west_east")

        def put(name, dims, data, dtype="f4"):
            if name in args.omit:
                return
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
        for name, value in (("PSFC", 100000.0), ("MUB", 10000.0), ("MAPFAC_M", 1.0),
                            ("SST", 290.0), ("TSK", 290.0), ("LANDMASK", 1.0)):
            put(name, xy, value)
        put("MAPFAC_U", ("Time", "south_north", "west_east_stag"), 1.0)
        put("MAPFAC_V", ("Time", "south_north_stag", "west_east"), 1.0)
        put("C1H", ("Time", "bottom_top"), 1.0)
        put("C2H", ("Time", "bottom_top"), 0.0)
        put("RDNW", ("Time", "bottom_top"), 1.0 / NZ)
        put("XLAT", xy, 30.0)
        put("XLONG", xy, -100.0)
        put("IVGTYP", xy, 16, "i4")
        put("ISLTYP", xy, 14, "i4")

        for name, data in surface_fields().items():
            for bad_name, value in args.fill:
                if bad_name == name:
                    data = np.full_like(data, float(value))
            for bad_name, value in args.bad:
                if bad_name == name:
                    data = data.copy()
                    # i = 8 is the first column of the deck's second box: the second rank
                    # of a two-rank run holds it, and it is a ghost column of the first box
                    data[..., 3, 8] = float(value)
            put(name, soil if data.ndim == 3 else xy, data[None, ...])
        zero = 0.0 if args.zero_layers else 1.0
        put("ZS", ("Time", "soil_layers_stag"), zero * (ZS_SURFACE if args.surface_zs else ZS)[None, :])
        put("DZS", ("Time", "soil_layers_stag"), zero * DZS[None, :])

    with nc.Dataset("wrflowinp_d01", "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("Time", 1)
        ds.createDimension("DateStrLen", len(DATE))
        ds.createDimension("south_north", NY)
        ds.createDimension("west_east", NX)
        ds.createVariable("Times", "S1", ("Time", "DateStrLen"))[:] = np.array(list(DATE), dtype="S1")
        ds.createVariable("SST", "f4", ("Time", "south_north", "west_east"))[:] = 290.0


if __name__ == "__main__":
    main()
