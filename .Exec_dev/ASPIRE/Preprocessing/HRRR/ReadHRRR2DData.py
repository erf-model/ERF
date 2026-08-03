import pygrib
import numpy as np
import struct
from pyproj import Proj, Transformer, CRS
import matplotlib.pyplot as plt
import sys
import os
from scipy.interpolate import interp1d
from collections import Counter
import subprocess
from pathlib import Path
from mpi4py import MPI

from IO import write_binary_vtk_structured_grid
from PlotImages import PlotVar

const_g = 9.81

def make_png_name(output_dir, var, date_time_forecast_str):

    output_png = output_dir + var + "_" + date_time_forecast_str + ".png"
    return output_png

def ReadHRRR2DData(file_path, area, lambert_conformal):
    # 1. Define master list of variable keys once
    VAR_KEYS = ["prate", "apcp", "sbcape", "mlcape", "mucape", "pwat", "ustm", "vstm", "cin", "st_hel","comp_rad_ref"]

    # 2. Dynamically build fields and units dictionaries
    fields = {key: None for key in VAR_KEYS}
    units = {key: "unknown" for key in VAR_KEYS}
    varnames = {key: "unknown" for key in VAR_KEYS}

    lats, lons = None, None
    printed_time = False
    date_time_forecast_str = ""

    with pygrib.open(file_path) as grbs:
        for grb in grbs:
            if not printed_time:
                year, month, day, hour = (
                    grb.year,
                    grb.month,
                    grb.day,
                    grb.hour,
                )
                minute = getattr(grb, "minute", 0)
                date_time_forecast_str = (
                    f"{year:04d}_{month:02d}_{day:02d}_{hour:02d}_{minute:02d}"
                )
                print(f"Datetime string: {date_time_forecast_str}")
                printed_time = True

            #print(f"Variable: {grb.name}, Units: {grb.parameterUnits}")
            # Extract 2D Fields
            if grb.name == "Total Precipitation":
                tmp = "apcp"
                varnames[tmp] = grb.name
                fields[tmp]   = grb.values
                units[tmp]    = grb.units
            elif grb.name == "Precipitation rate":
                tmp = "prate"
                varnames[tmp] = grb.name
                fields[tmp]  = grb.values
                units[tmp]    = grb.units
            elif grb.name == "Convective available potential energy":
                if grb.level == 0:
                    tmp = "sbcape"
                    varnames[tmp] = grb.name
                    fields[tmp]   = grb.values
                    units[tmp]    = grb.units
                elif grb.level == 9000:
                    tmp = "mlcape"
                    varnames[tmp] = grb.name
                    fields[tmp]   = grb.values
                    units[tmp]    = grb.units
                elif grb.level == 25500:
                    tmp = "mucape"
                    varnames[tmp] = grb.name
                    fields[tmp]   = grb.values
                    units[tmp]    = grb.units
            elif grb.name == "Precipitable water":
                tmp = "pwat"
                varnames[tmp] = grb.name
                fields[tmp]   = grb.values
                units[tmp]    = grb.units
            elif grb.name == "U-component storm motion":
                tmp = "ustm"
                varnames[tmp] = grb.name
                fields[tmp]   = grb.values
                units[tmp]    = grb.units
            elif grb.name == "V-component storm motion":
                tmp = "vstm"
                varnames[tmp] = grb.name
                fields[tmp]   = grb.values
                units[tmp]    = grb.units
            elif grb.name == "Convective inhibition":
                tmp = "cin"
                varnames[tmp] = grb.name
                fields[tmp]   = grb.values
                units[tmp]    = grb.units
            elif grb.name == "Storm relative helicity":
                tmp = "st_hel"
                varnames[tmp] = grb.name
                fields[tmp]   = grb.values
                units[tmp]    = grb.units
            elif grb.name == "Maximum/Composite radar reflectivity":
                tmp = "comp_rad_ref"
                varnames[tmp] = grb.name
                fields[tmp]   = grb.values
                units[tmp]    = grb.units

            if lats is None or lons is None:
                lats, lons = grb.latlons()

    # Bounding Box: area = [lat_max, lon_min, lat_min, lon_max]
    lat_max, lon_min, lat_min, lon_max = area[0], area[1], area[2], area[3]

    # Transform native 2D grids to projection coordinates
    transformer = Transformer.from_crs(
        "EPSG:4326", lambert_conformal, always_xy=True
    )
    x_native, y_native = transformer.transform(lons, lats)

    # Project domain limits
    x_min_b, y_min_b = transformer.transform(lon_min, lat_min)
    x_max_b, y_max_b = transformer.transform(lon_max, lat_max)

    # Mask native grid within the targeted area
    mask = (
        (x_native >= min(x_min_b, x_max_b))
        & (x_native <= max(x_min_b, x_max_b))
        & (y_native >= min(y_min_b, y_max_b))
        & (y_native <= max(y_min_b, y_max_b))
    )

    i_indices, j_indices = np.where(mask)
    lat_start, lat_end = i_indices.min(), i_indices.max()
    lon_start, lon_end = j_indices.min(), j_indices.max()

    # Slice spatial coordinates and variables consistently
    x_grid = x_native[lat_start : lat_end + 1, lon_start : lon_end + 1]
    y_grid = y_native[lat_start : lat_end + 1, lon_start : lon_end + 1]

    output_fields = {}
    for key in fields:
        if fields[key] is not None:
            output_fields[key] = fields[key][
                lat_start : lat_end + 1, lon_start : lon_end + 1
            ]


    output_fields["prate"] =  output_fields["prate"]*3600.0
    units["prate"] = "mm/hr"

    # Plot output loop over all variables
    for var_dir_name in VAR_KEYS:
        output_dir = os.path.join("./Output/Thunderstorm", var_dir_name)
        os.makedirs(output_dir, exist_ok=True)

        output_png = make_png_name(
            output_dir + "/", var_dir_name, date_time_forecast_str
        )

        PlotVar(
            varnames[var_dir_name],
            output_fields[var_dir_name],
            units[var_dir_name],
            date_time_forecast_str,
            x_grid,
            y_grid,
            lambert_conformal,
            output_png,
        )
