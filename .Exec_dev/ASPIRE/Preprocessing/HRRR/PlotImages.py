import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np


def PlotMaxDBZ_LCC(max_dbz, x_grid, y_grid, lambert_proj, output_file):

    # Cartopy CRS from your pyproj string
    lcc = ccrs.LambertConformal(
        central_longitude=float(
            lambert_proj.split("+lon_0=")[1].split()[0]
        ),
        central_latitude=float(
            lambert_proj.split("+lat_0=")[1].split()[0]
        ),
        standard_parallels=(
            float(lambert_proj.split("+lat_1=")[1].split()[0]),
            float(lambert_proj.split("+lat_2=")[1].split()[0])
        )
    )

    fig = plt.figure(figsize=(10, 8))

    ax = plt.axes(projection=lcc)

    # Domain
    ax.set_extent([
        x_grid.min(),
        x_grid.max(),
        y_grid.min(),
        y_grid.max()
    ], crs=lcc)


    # Plot reflectivity
    levels = np.arange(0, 75, 5)

    cf = ax.contourf(
        x_grid,
        y_grid,
        max_dbz,
        levels=levels,
        cmap="turbo",
        extend="max",
        transform=lcc
    )


    # State boundaries
    ax.add_feature(
        cfeature.STATES,
        linewidth=0.5
    )

    ax.add_feature(
        cfeature.COASTLINE,
        linewidth=0.7
    )

    ax.add_feature(
        cfeature.BORDERS,
        linewidth=0.5
    )


    # Color bar
    cbar = plt.colorbar(
        cf,
        ax=ax,
        shrink=0.8
    )

    cbar.set_label(
        "Maximum Simulated Reflectivity (dBZ)"
    )


    plt.title(
        "HRRR Maximum Simulated Reflectivity"
    )

    plt.savefig(
        output_file,
        dpi=200,
        bbox_inches="tight"
    )

    plt.close()
