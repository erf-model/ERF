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
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="lightgray")
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white")
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.8)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.8, edgecolor="black")
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.5, edgecolor="black")

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
    track = np.loadtxt(track_file)
    if track.ndim != 2 or track.shape[1] < 2:
        raise ValueError("Track file must have at least two columns (lon, lat)")
    lons, lats = track[:, 0], track[:, 1]
    ax.plot(lons, lats, "-o", color=color, linewidth=2, markersize=1,
            transform=ccrs.Geodetic(), label=label)
    return ax


def plot_actual_hurricane_track(ax, actual_track_file, color="black",
                                label="Actual track", marker="-o", linewidth=2):
    track = np.loadtxt(actual_track_file)
    if track.ndim != 2 or track.shape[1] < 2:
        raise ValueError("Track file must have at least two columns (lon, lat)")
    lons, lats = track[:, 0], track[:, 1]
    ax.plot(lons, lats, marker, color=color, linewidth=linewidth, markersize=1,
            transform=ccrs.Geodetic(), label=label)
    return ax


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot lat/lon map with optional hurricane track and diagnostics."
    )

    parser.add_argument("--area", required=True,
                        help="Comma-separated lon_min,lon_max,lat_min,lat_max")
    parser.add_argument("--erf_track", default=None)
    parser.add_argument("--actual_track", default=None)
    parser.add_argument("--outfile_track", default="map.png")

    # --- Max velocity ---
    parser.add_argument("--erf_maxvel", default=None)
    parser.add_argument("--actual_maxvel", default=None)
    parser.add_argument("--wrf_maxvel", default=None)
    parser.add_argument("--outfile_maxvel", default=None)

    # --- Min pressure ---
    parser.add_argument("--erf_minpressure", default=None)
    parser.add_argument("--actual_minpressure", default=None)
    parser.add_argument("--wrf_minpressure", default=None)
    parser.add_argument("--outfile_minpressure", default=None)

    args = parser.parse_args()

    # ---------------- Map plot ----------------
    lon_min, lon_max, lat_min, lat_max = map(float, args.area.split(","))

    fig, ax = plt.subplots(figsize=(10, 8),
                           subplot_kw={"projection": ccrs.PlateCarree()})
    plot_latlon_map(ax, lon_min, lon_max, lat_min, lat_max)

    if args.erf_track:
        plot_hurricane_track(ax, args.erf_track)
    if args.actual_track:
        plot_actual_hurricane_track(ax, args.actual_track)

    ax.legend(loc="upper right")
    plt.title("Latitude/Longitude Map", fontsize=16)
    plt.savefig(args.outfile_track, dpi=150, bbox_inches="tight")
    plt.close()

    # ---------------- Max velocity plot ----------------
    maxvel_files = [args.erf_maxvel, args.actual_maxvel, args.wrf_maxvel]
    if any(f is not None for f in maxvel_files):
        if args.outfile_maxvel is None:
            raise ValueError("--outfile_maxvel must be provided if any maxvel input is given")

        plt.figure(figsize=(8, 5))

        if args.erf_maxvel:
            data = np.loadtxt(args.erf_maxvel)
            data_s = savgol_filter(data[:, 1], window_length=24, polyorder=3)
            plt.plot(data[:, 0], data_s / 1.852, "-xr", linewidth=2, label="ERF")

        if args.wrf_maxvel:
            data = np.loadtxt(args.wrf_maxvel)
            plt.plot(data[:, 0], data[:, 1], "-ob", linewidth=2, label="WRF")

        if args.actual_maxvel:
            data = np.loadtxt(args.actual_maxvel)
            plt.plot(data[:, 0], data[:, 1], "-k", linewidth=2, label="Actual")

        plt.xlabel("Time (hours)")
        plt.ylabel("Max wind speed (km/hr)")
        plt.legend()
        plt.savefig(args.outfile_maxvel, dpi=300, bbox_inches="tight")
        plt.close()

    # ---------------- Min pressure plot ----------------
    minp_files = [args.erf_minpressure,
                  args.actual_minpressure,
                  args.wrf_minpressure]

    if any(f is not None for f in minp_files):
        if args.outfile_minpressure is None:
            raise ValueError("--outfile_minpressure must be provided if any minpressure input is given")

        plt.figure(figsize=(8, 5))

        if args.erf_minpressure:
            data = np.loadtxt(args.erf_minpressure)
            plt.plot(data[:, 0], data[:, 1], "-xr", linewidth=2, label="ERF")

        if args.wrf_minpressure:
            data = np.loadtxt(args.wrf_minpressure)
            plt.plot(data[:, 0], data[:, 1], "-ob", linewidth=2, label="WRF")

        if args.actual_minpressure:
            data = np.loadtxt(args.actual_minpressure)
            plt.plot(data[:, 0], data[:, 1], "-k", linewidth=2, label="Actual")

        plt.xlabel("Time (hours)")
        plt.ylabel("Minimum pressure")
        plt.legend()
        plt.savefig(args.outfile_minpressure, dpi=300, bbox_inches="tight")
        plt.close()


