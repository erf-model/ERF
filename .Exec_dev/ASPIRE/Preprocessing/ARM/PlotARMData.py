#!/usr/bin/env python3

import argparse
import glob
import os

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature


# ------------------------------------------------------
# Command line arguments
# ------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Plot ARM precipitation data and station location"
)

parser.add_argument(
    "--arm-data-dir",
    required=True,
    help="Directory containing ARM .cdf files"
)

args = parser.parse_args()

arm_data_dir = args.arm_data_dir


# ------------------------------------------------------
# Find ARM files
# ------------------------------------------------------

files = sorted(
    glob.glob(os.path.join(arm_data_dir, "*.cdf"))
)

if len(files) == 0:
    raise RuntimeError(
        f"No ARM .cdf files found in {arm_data_dir}"
    )

print("Reading ARM files:")
for f in files:
    print("  ", f)


# ------------------------------------------------------
# Read ARM data
# ------------------------------------------------------

ds = xr.open_mfdataset(
    files,
    combine="by_coords"
)

for var in ds.data_vars:
    long_name = ds[var].attrs.get("long_name", "No long name")
    units = ds[var].attrs.get("units", "")
    print(f"{var}: {long_name} ({units})")

print(ds)


# ------------------------------------------------------
# Extract station location
# ------------------------------------------------------

lat = float(ds["lat"].values.flat[0])
lon = float(ds["lon"].values.flat[0])
alt = float(ds["alt"].values.flat[0])

print("\nARM Station:")
print(f" Latitude  : {lat:.5f}")
print(f" Longitude : {lon:.5f}")
print(f" Altitude  : {alt:.1f} m")

def plot_function(varname, units):
    plt.figure(figsize=(12, 5))
    ds[varname].plot()
    plt.title(
    "ARM Tipping Bucket Rain Gauge\n"
    f"{str(ds.time.values[0])[:10]} to "
    f"{str(ds.time.values[-1])[:10]}"
    )
    plt.xlabel("Time")
    ylabel_str = varname  + "(" + units +")"
    plt.ylabel(ylabel_str)
    plt.grid(True)
    plt.tight_layout()
    plt_filename = varname + ".png"
    plt.savefig(
    plt_filename,
    dpi=300,
    bbox_inches="tight"
    )
    print(f"Saved {plt_filename:s}")    

# ------------------------------------------------------
# Plot accumulated precipitation
# ------------------------------------------------------

var = 'tbrg_precip_total_corr';
# Cumulative sum to get running total across the period
ds['precip_accum_cumsum'] = ds['tbrg_precip_total_corr'].cumsum()
ds['precip_accum_cumsum'].plot()
units = ds[var].attrs.get("units", "")
plot_function("precip_accum_cumsum", units)

# ------------------------------------------------------
# Plot wind speed arithmetic mean
# ------------------------------------------------------
var = 'wspd_arith_mean';
units = ds[var].attrs.get("units", "");
plot_function(var, units)

# ------------------------------------------------------
# Plot temperature mean
# ------------------------------------------------------
var = 'temp_mean';
units = ds[var].attrs.get("units", "");
plot_function(var, units)

# ------------------------------------------------------
# Plot temperature mean
# ------------------------------------------------------
var = 'rh_mean';
units = ds[var].attrs.get("units", "");
plot_function(var, units)


# ------------------------------------------------------
# Plot station location (~300 km x 300 km region)
# ------------------------------------------------------

fig = plt.figure(figsize=(7, 7))

ax = plt.axes(
    projection=ccrs.PlateCarree()
)

# 300 km x 300 km box centered on station
half_box_km = 1000.0

lat_extent = half_box_km / 111.0

lon_extent = (
    half_box_km /
    (111.0 * np.cos(np.deg2rad(lat)))
)

ax.set_extent(
    [
        lon - lon_extent,
        lon + lon_extent,
        lat - lat_extent,
        lat + lat_extent
    ],
    crs=ccrs.PlateCarree()
)


# Map features

ax.add_feature(
    cfeature.STATES,
    linewidth=0.5
)

ax.add_feature(
    cfeature.BORDERS,
    linewidth=0.8
)

ax.add_feature(
    cfeature.COASTLINE,
    linewidth=0.8
)


# Station marker

ax.plot(
    lon,
    lat,
    marker="*",
    markersize=18,
    transform=ccrs.PlateCarree()
)

ax.text(
    lon + 0.05,
    lat + 0.05,
    "ARM Station",
    transform=ccrs.PlateCarree(),
    fontsize=12
)


plt.title(
    "ARM Station Location\n"
    f"Lat={lat:.3f}, Lon={lon:.3f}"
)

plt.savefig(
    "ARM_station_location.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Saved ARM_station_location.png")


