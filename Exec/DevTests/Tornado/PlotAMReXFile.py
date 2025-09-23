import yt
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import argparse

def slice_and_contour_latlon(plotfile, field, axis="z", coord=0.0,
                             outfile_data="slice.npy", outfile_plot="contour.png",
                             nx=512, ny=512,
                             lon_min=-100, lon_max=-80,
                             lat_min=30, lat_max=46,
                             loc_lat=None, loc_lon=None, loc_name=None):
    """
    Extract a 2D slice from an AMReX plotfile and plot in lat/lon coordinates with Cartopy.
    """

    # Fixed scaling factor for latitude axis
    y_scale = 1.05

    # Load dataset
    ds = yt.load(plotfile)

    # Make slice
    data_source = ds.slice(axis, coord)

    # Get FRB (Fixed Resolution Buffer)
    frb = data_source.to_frb((1.0, "unitary"), (nx, ny))
    arr = np.array(frb[field])

    # Save raw data
    np.save(outfile_data, arr)
    print(f"Slice data saved to {outfile_data}")

    # Normalized coordinates
    x_unit = np.linspace(float(frb.bounds[0]), float(frb.bounds[1]), nx)
    y_unit = np.linspace(float(frb.bounds[2]), float(frb.bounds[3]), ny)

    # Map to lon/lat with adjustable y scaling
    lon = lon_min + (lon_max - lon_min) * (x_unit - x_unit.min()) / (x_unit.max() - x_unit.min())
    lat = lat_min + (lat_max - lat_min) * (y_unit - y_unit.min()) / (y_unit.max() - y_unit.min())
    lat = lat_min + (lat - lat_min) * y_scale  # apply scaling

    LON, LAT = np.meshgrid(lon, lat)

    # Set up Cartopy plot
    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_min + (lat_max - lat_min) * y_scale],
                  crs=ccrs.PlateCarree())

    # Filled contour
    vmin, vmax = 5.0, 50.0
    cs = ax.contourf(LON, LAT, arr, levels=np.arange(vmin, vmax+5, 5),
                     cmap="rainbow", extend="max", transform=ccrs.PlateCarree())
    cbar = plt.colorbar(cs, ax=ax, orientation="vertical", shrink=0.95)
    cbar.set_label(field)

    # Add map features
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor='black')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.8, edgecolor='black')
    ax.coastlines('50m', linewidth=0.8)

    # Add gridlines and lat/lon ticks
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      x_inline=False, y_inline=False, linewidth=0.33, color='k', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"fontsize": 15}
    gl.ylabel_style = {"fontsize": 15}

    # Add location marker if provided
    if loc_lat is not None and loc_lon is not None:
        ax.plot(loc_lon, loc_lat, color='black', marker='*', markersize=8,
                transform=ccrs.Geodetic())
        if loc_name is not None:
            ax.text(loc_lon + 0.3, loc_lat + 0.3, loc_name,
                    fontsize=8, color='black', ha='left', transform=ccrs.Geodetic())

    plt.axis("tight")
    ax.set_ylim(top=lat_max)
    plt.savefig(outfile_plot, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Contour plot saved to {outfile_plot}")

    return arr


# ========================
# Command-line interface
# ========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot ERF slice in lat/lon")
    parser.add_argument("plotfile", help="Path to AMReX plotfile")
    parser.add_argument("--var", required=True, help="Variable name to plot")
    parser.add_argument("--area", required=True,
                        help="Comma-separated lon_min,lon_max,lat_min,lat_max (e.g. --area=-103.5,-80.5,29.8,46.18)")
    args = parser.parse_args()

    # Parse area
    lon_min, lon_max, lat_min, lat_max = map(float, args.area.split(","))

    # Filenames based on var
    outfile_data = f"{args.var}_slice.npy"
    outfile_plot = f"{args.var}_slice.png"

    slice_and_contour_latlon(
        plotfile=args.plotfile,
        field=args.var,
        axis="z",
        coord=0.0,
        outfile_data=outfile_data,
        outfile_plot=outfile_plot,
        nx=512, ny=512,
        lon_min=lon_min, lon_max=lon_max,
        lat_min=lat_min, lat_max=lat_max,
        loc_lat=37.0838, loc_lon=-94.5132, loc_name="Joplin"
    )

