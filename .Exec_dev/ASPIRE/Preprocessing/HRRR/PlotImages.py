import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

def PlotVar(
    varname,
    vardata,
    units,
    date_time_forecast_str,
    x_grid,
    y_grid,
    lambert_proj,
    output_file,
):
    # 1. Cartopy CRS from pyproj string
    lcc = ccrs.LambertConformal(
        central_longitude=float(
            lambert_proj.split("+lon_0=")[1].split()[0]
        ),
        central_latitude=float(lambert_proj.split("+lat_0=")[1].split()[0]),
        standard_parallels=(
            float(lambert_proj.split("+lat_1=")[1].split()[0]),
            float(lambert_proj.split("+lat_2=")[1].split()[0]),
        ),
    )

    # Use explicit subplots to tie axes cleanly to this specific figure instance
    fig, ax = plt.subplots(
        figsize=(10, 8), subplot_kw={"projection": lcc}
    )

    # Domain extent
    ax.set_extent(
        [x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()], crs=lcc
    )

    # --- Robust Level Generation ---
    vmin = np.nanmin(vardata)
    vmax = np.nanmax(vardata)

    # Handle NaNs (entire array is NaN or masked) or uniform zero/flat fields
    if np.isnan(vmin) or np.isnan(vmax) or vmin == vmax:
        if np.isnan(vmin):
            vmin, vmax = 0.0, 1.0
        else:
            vmin = vmin - 0.1
            vmax = vmax + 0.1

    levels = np.linspace(vmin, vmax, 15)

    # Guarantee strictly increasing, unique level values
    if not np.all(np.diff(levels) > 0):
        levels = np.unique(levels)
        if len(levels) < 2:
            levels = np.array([vmin - 0.1, vmax + 0.1])
    # -------------------------------

    cf = ax.contourf(
        x_grid,
        y_grid,
        vardata,
        levels=levels,
        cmap="turbo",
        extend="both",  # 'both' handles under/over values safely
        transform=lcc,
    )

    # Add geographic features
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.7)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.5)

    # Colorbar bound strictly to the local figure (`fig.colorbar`)
    cbar = fig.colorbar(cf, ax=ax, shrink=0.8)
    cbar.set_label(f"{varname} [{units}]", fontsize=10, labelpad=10)

    # Bind title explicitly to local axes (prevents global title overlaps)
    title_text = f"{varname} {date_time_forecast_str}"
    ax.set_title(title_text, fontsize=12)

    # Save figure
    fig.savefig(output_file, dpi=200, bbox_inches="tight")

    # --- COMPLETE RESET & CLEANUP ---
    plt.close(fig)  # Closes ALL active figure windows globally
