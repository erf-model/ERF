import pygrib
import numpy as np
import struct
from pyproj import Proj, Transformer, CRS
import matplotlib.pyplot as plt
import sys
import os
from scipy.interpolate import interp1d
from collections import Counter

from IO import write_binary_vtk_structured_grid

from Thunderstorm_Squall_HeavyRain import compute_Thunderstorm_Squall_HeavyRain
from PlotImages import PlotMaxDBZ_LCC

const_g = 9.81

def remove_extra_levels(data, levels, target_levels):

    levels = np.array(levels)
    target_levels = np.array(target_levels)

    indices = []

    for lev in target_levels:
        idx = np.where(levels == lev)[0]
        if len(idx) > 0:
            indices.append(idx[0])

    indices = np.array(indices)

    return data[indices,:,:], levels[indices]

def filter_fields(fields, levels, pressure_levels):

    for key in fields:
        fields[key], levels[key] = remove_extra_levels(
            fields[key],
            levels[key],
            pressure_levels
        )

    return fields, levels

def ReadHRRR3DData(file_path, area, lambert_conformal):

    fields = {
    "temp": [],
    "ght": [],
    "uvel": [],
    "vvel": [],
    "wvel": [],
    "qv": [],
    "qc": [],
    "qr": [],
    "qs": [],
    "qg": [],
    "rh": [],
    }
    theta_3d_hr3 = []
    pressure_3d_hr3 = []
    lats, lons = None, None  # To store latitude and longitude grids

    levels = {
    "temp": [],
    "ght": [],
    "uvel": [],
    "vvel": [],
    "wvel": [],
    "qv": [],
    "qc": [],
    "qr": [],
    "qs": [],
    "qg": [],
    "rh": [],
    }

    printed_time = False

    date_time_forecast_str = ""

    #with pygrib.open(file_path) as grbs:
    #    for grb in grbs:
    #        if "Geopotential" in grb.name:
    #            levels.append(grb.level)

    #print(sorted(set(levels)))
    #print(Counter(levels))
    #exit()

    with pygrib.open(file_path) as grbs:
        for grb in grbs:

            if not printed_time:
                year = grb.year
                month = grb.month
                day = grb.day
                hour = grb.hour

                minute = grb.minute if hasattr(grb, 'minute') else 0
                print(f"Date: {year}-{month:02d}-{day:02d}, Time: {hour:02d}:{minute:02d} UTC")
                date_time_forecast_str = f"{year:04d}_{month:02d}_{day:02d}_{hour:02d}_{minute:02d}"
                print(f"Datetime string: {date_time_forecast_str}")
                printed_time = True

            #print(f"Variable: {grb.name}, Level: {grb.level}, Units: {grb.parameterUnits}")
            if "Temperature" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                # Append temperature values
                fields["temp"].append(grb.values)
                levels["temp"].append(grb.level)

            if "Geopotential" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["ght"].append(grb.values)
                levels["ght"].append(grb.level)

            if "U component of wind" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["uvel"].append(grb.values)
                levels["uvel"].append(grb.level)

            if "V component of wind" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["vvel"].append(grb.values)
                levels["vvel"].append(grb.level)

            if "Vertical velocity" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["wvel"].append(grb.values)
                levels["wvel"].append(grb.level)

            if "Specific humidity" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["qv"].append(grb.values)
                levels["qv"].append(grb.level)

            if "Cloud mixing ratio" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["qc"].append(grb.values)
                levels["qc"].append(grb.level)

            if "Rain mixing ratio" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["qr"].append(grb.values)
                levels["qr"].append(grb.level)

            if "Snow mixing ratio" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["qs"].append(grb.values)
                levels["qs"].append(grb.level)

            if "Graupel (snow pellets)" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["qg"].append(grb.values)
                levels["qg"].append(grb.level)

            if "Relative humidity" in grb.name and grb.typeOfLevel == "isobaricInhPa":
                fields["rh"].append(grb.values)
                levels["rh"].append(grb.level)

            # Retrieve latitude and longitude grids (once)
            if lats is None or lons is None:
                lats, lons = grb.latlons()

    for key in fields:
        if len(fields[key]) > 0:
            fields[key] = np.stack(fields[key], axis=0)
            levels[key] = np.array(levels[key])
        else:
            print(f"Warning: {key} not found in file")

    #pressure_3d_hr3 = np.stack(pressure_3d_hr3, axis=0)
    # Get the size of each dimension
    dim1, dim2, dim3 = fields["uvel"].shape

    # Print the sizes
    print(f"Size of dimension 1: {dim1}")
    print(f"Size of dimension 2: {dim2}")
    print(f"Size of dimension 3: {dim3}")


    # Extract unique latitude and longitude values
    unique_lats = np.unique(lats[:, 0])  # Take the first column for unique latitudes
    unique_lons = np.unique(lons[0, :])  # Take the first row for unique longitudes

    print("Min max lat lons are ", unique_lats[0], unique_lats[-1], unique_lons[0], unique_lons[-1]);

    nlats = len(unique_lats)
    nlons = len(unique_lons)

    lat_max = area[0]
    lon_min = area[1]
    lat_min = area[2]
    lon_max = area[3]

    print("Lat/lon min/max are ", lat_min, lat_max, lon_min, lon_max, flush=True)

    # Example: regular grid
    lat_resolution = unique_lats[1] - unique_lats[0]
    lon_resolution = unique_lons[1] - unique_lons[0]

    lat_start = int((lat_min - unique_lats[0]) / lat_resolution)
    lat_end   = int((lat_max - unique_lats[0]) / lat_resolution)
    lon_start = int((lon_min - unique_lons[0]) / lon_resolution)
    lon_end   = int((lon_max - unique_lons[0]) / lon_resolution)

    domain_lats = unique_lats[lat_start:lat_end+1]
    domain_lons = unique_lons[lon_start:lon_end+1]

    print("The min max are",(lat_start, lat_end, lon_start, lon_end));

    nx = domain_lons.shape[0]
    ny = domain_lats.shape[0]

    print("nx and ny here are ", nx, ny)

    pressure_levels = levels["qc"]

    fields, levels = filter_fields(fields, levels, pressure_levels)

    for key in fields:
        print(
        key,
        fields[key].shape,
        levels[key]
        )

    nz = len(pressure_levels)

    output_fields = {}

    for key in ["ght", "uvel", "vvel", "wvel",
            "qv", "qc", "qr", "qs", "qg", "rh", "temp"]:
        output_fields[key] = np.zeros((ny, nx, nz))

    # Create meshgrid
    x_grid, y_grid = np.meshgrid(domain_lons, domain_lats)
    lon_grid, lat_grid = np.meshgrid(domain_lons, domain_lats)

    #lambert_conformal = CRS.from_proj4(
    #"+proj=lcc +lat_1=30 +lat_2=60 +lat_0=38.5 +lon_0=-97 +datum=WGS84 +units=m +no_defs")

    transformer = Transformer.from_crs("EPSG:4326", lambert_conformal, always_xy=True)

    # Convert the entire grid to UTM
    x_grid, y_grid = transformer.transform(lon_grid, lat_grid)

    for k in range(nz):
        for key in output_fields:
            output_fields[key][:,:,k] = fields[key][k,lat_start:lat_end+1,lon_start:lon_end+1]

    z_grid   = output_fields["ght"]

    for k, p in enumerate(levels["ght"]):
        ght = fields["ght"][k]

        print(
            f"{p:4d} hPa : "
            f"min={ght.min():10.2f} "
            f"max={ght.max():10.2f}"
        )

    Rd = 287.05  # J/(kg K)

    rho_air = np.zeros_like(output_fields["temp"])

    for k, p_hpa in enumerate(levels["temp"]):
        p_pa = p_hpa * 100.0
        rho_air[:, :, k] = p_pa / (Rd * output_fields["temp"][:, :, k])

    # Compute event specific quantities
    max_dbz_2d, max_dbz_3d = compute_Thunderstorm_Squall_HeavyRain (nz, rho_air,
                                       output_fields["temp"],
                                       output_fields["qr"],
                                       output_fields["qs"],
                                       output_fields["qg"],)

    #print(output_fields["ght"][:,:,39])

    scalars = {
         "uvel": output_fields["uvel"],
         "vvel": output_fields["vvel"],
         "wvel": output_fields["wvel"],
         "qv": output_fields["qv"],
         "qc": output_fields["qc"],
         "qr": output_fields["qr"],
         "qs": output_fields["qs"],
         "qg": output_fields["qg"],
         "rh": output_fields["rh"],
         "temp": output_fields["temp"],
         "max_dbz": max_dbz_3d,
    }

    dir_path = "Output"
    os.makedirs(dir_path, exist_ok=True)

     # Extract the filename from the full path
    filename = os.path.basename(file_path)

    output_vtk = "./Output/VTK/3D/HRRRDomain/HRRR_" + date_time_forecast_str + ".vtk"

    #write_binary_vtk_structured_grid(output_vtk, x_grid, y_grid, z_grid, scalars)

    output_png = "./Output/Thunderstorm/max_dbz_" + date_time_forecast_str + ".png"

    PlotMaxDBZ_LCC(
    max_dbz_2d,
    x_grid,
    y_grid,
    lambert_conformal,
    output_png)

