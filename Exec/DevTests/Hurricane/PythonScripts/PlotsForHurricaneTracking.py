#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import argparse
from scipy.signal import savgol_filter


# ------------------------------
# Function 1: plot the base map
# ------------------------------
def plot_latlon_map(ax, lon_min, lon_max, lat_min, lat_max):
    """
    Draw a WRF-style map with coastlines, borders, and states on a Cartopy axis.
    """
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Map features
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="lightgray")
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white")
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.8)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.8, edgecolor="black")
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.5, edgecolor="black")

    # Gridlines & labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      x_inline=False, y_inline=False,
                      linewidth=0.33, color="k", alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"fontsize": 13}
    gl.ylabel_style = {"fontsize": 13}

    return ax


# --------------------------------
# Function 2: plot hurricane track
# --------------------------------
def plot_hurricane_track(ax, track_file, color="red", label="ERF track"):
    """
    Plot a hurricane track from a text file containing lon/lat pairs on an existing Cartopy axis.
    """
    track = np.loadtxt(track_file)
    if track.ndim != 2 or track.shape[1] < 2:
        raise ValueError("Track file must have at least two columns (lon, lat)")
    lons, lats = track[:, 0], track[:, 1]
    ax.plot(lons, lats, "-o", color=color, linewidth=2, markersize=1,
            transform=ccrs.Geodetic(), label=label)
    return ax


def plot_actual_hurricane_track(ax, actual_track_file, color="black", label="Actual track", marker="-o", linewidth=2):
    """
    Plot a hurricane track from a text file containing lon/lat pairs on an existing Cartopy axis.
    """
    track = np.loadtxt(actual_track_file)
    if track.ndim != 2 or track.shape[1] < 2:
        raise ValueError("Track file must have at least two columns (lon, lat)")
    lons, lats = track[:, 0], track[:, 1]
    ax.plot(lons, lats, marker, color=color, linewidth=linewidth, markersize=1,
            transform=ccrs.Geodetic(), label=label)
    return ax

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot lat/lon map with optional hurricane track and max wind speed overlay."
    )
    parser.add_argument("--area", required=True,
                        help="Comma-separated lon_min,lon_max,lat_min,lat_max "
                             "(e.g. --area=-103.5,-80.5,29.8,46.18)")
    parser.add_argument("--erf_track", default=None,
                        help="Path to text file with ERF hurricane track (lon lat columns)")
    parser.add_argument("--actual_track", default=None,
                        help="Path to text file with Actual hurricane track (lon lat columns)")
    parser.add_argument("--outfile_track", default="map.png",
                        help="Output image filename for the map")
    
    # Optional max velocity inputs
    parser.add_argument("--erf_maxvel", default=None, help="ERF maxvel vs time data")
    parser.add_argument("--actual_maxvel", default=None, help="Actual maxvel vs time data")
    parser.add_argument("--wrf_maxvel", default=None, help="WRF maxvel vs time data")
    parser.add_argument("--outfile_maxvel", default=None, help="Output image filename for maxvel plot")

    args = parser.parse_args()

    # --- Parse map area ---
    lon_min, lon_max, lat_min, lat_max = map(float, args.area.split(","))

    # --- Plot lat/lon map ---
    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": ccrs.PlateCarree()})
    plot_latlon_map(ax, lon_min, lon_max, lat_min, lat_max)

    # Plot hurricane tracks if provided
    if args.erf_track:
        plot_hurricane_track(ax, args.erf_track)
    if args.actual_track:
        plot_actual_hurricane_track(ax, args.actual_track)

    ax.legend(loc="upper right")
    plt.title("Latitude/Longitude Map", fontsize=16)
    plt.savefig(args.outfile_track, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure saved to {args.outfile_track}")

    # --- Max velocity plots ---
    maxvel_files = [args.erf_maxvel, args.actual_maxvel, args.wrf_maxvel]
    any_maxvel = any(f is not None for f in maxvel_files)

    if any_maxvel:
        if args.outfile_maxvel is None:
            raise ValueError("If any of --erf_maxvel, --actual_maxvel, or --wrf_maxvel is provided, "
                             "--outfile_maxvel must also be provided.")

        plt.figure(figsize=(8, 5))

        if args.erf_maxvel:
            data_erf = np.loadtxt(args.erf_maxvel)
            data_erf_smooth = savgol_filter(data_erf[:, 1], window_length=24, polyorder=3)
            plt.plot(data_erf[:, 0], data_erf_smooth / 1.852, '-xr', linewidth=2, label="ERF")

        if args.wrf_maxvel:
            data_wrf = np.loadtxt(args.wrf_maxvel)
            plt.plot(data_wrf[:, 0], data_wrf[:, 1], '-ob', linewidth=2, label="WRF")

        if args.actual_maxvel:
            data_actual = np.loadtxt(args.actual_maxvel)
            plt.plot(data_actual[:, 0], data_actual[:, 1], '-k', linewidth=2, label="Actual")

        plt.xlabel("Time (hours)")
        plt.ylabel("Max. wind speed (km/hr)")
        plt.legend()
        plt.savefig(args.outfile_maxvel, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Figure saved to {args.outfile_maxvel}")


