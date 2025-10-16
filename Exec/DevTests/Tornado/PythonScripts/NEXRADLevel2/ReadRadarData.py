import pyart
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Geod
import glob
import os
import argparse

# --- Parse command-line arguments ---
parser = argparse.ArgumentParser(description="Plot NEXRAD radar PPI with a city marker")
parser.add_argument("--lat", type=float, required=False, help="City latitude")
parser.add_argument("--lon", type=float, required=False, help="City longitude")
parser.add_argument("--data_folder", type=str, required=True, help="Folder with NEXRAD .gz files")
args = parser.parse_args()

lat_city = args.lat
lon_city = args.lon
# Folder containing NEXRAD Level II files
data_folder = args.data_folder

# Create output folder
output_folder = os.path.join(data_folder, "Images")
os.makedirs(output_folder, exist_ok=True)

geod = Geod(ellps="WGS84")

# Get all *.gz files sorted
files = sorted(glob.glob(os.path.join(data_folder, "*.gz")))

for fpath in files:
    # Extract timestamp from filename
    fname = os.path.basename(fpath)
    # Example: KSGF20110522_223858_V03.gz -> timestamp = 20110522_223858
    timestamp = fname[4:19]

    print(f"Processing {fname} ...")


    # --- Read radar ---
    radar = pyart.io.read_nexrad_archive(fpath)

    for i in range(radar.nsweeps):
        sl = radar.get_slice(i)
        v = radar.fields['velocity']['data'][sl]
        n_valid = (~v.mask).sum() if hasattr(v, "mask") else v.size
        elev = radar.elevation['data'][sl][0]
        print(f"Sweep {i}: elevation={elev:.2f}°, valid gates={n_valid}")


    # --- Sweep selection ---
    sweep = 1
    slc = radar.get_slice(sweep)
    vel = radar.fields['velocity']['data'][slc]

    # Radar gate Cartesian coordinates
    x, y, z = radar.get_gate_x_y_z(sweep)  # meters
    x_flat = x.flatten()
    y_flat = y.flatten()
    z_flat = z.flatten()
    vel_flat = vel.flatten()
    mask = ~vel_flat.mask if hasattr(vel_flat, "mask") else np.ones_like(vel_flat, dtype=bool)
    x_flat = x_flat[mask]
    y_flat = y_flat[mask]
    z_flat = z_flat[mask]
    vel_flat = vel_flat[mask]

    # Convert to cylindrical coordinates
    r = np.hypot(x_flat, y_flat)       # radius (m)
    theta = np.arctan2(y_flat, x_flat) # azimuth angle (rad)
    phi = np.arctan2(z_flat, r)        # elevation angle (rad)

    # Print some sample values
    print("Sample cylindrical coordinates and velocity (radius m, theta deg, phi deg, vel m/s):")
    for i in range(0, min(10, len(r))):
        print(f"r={r[i]:.1f}, theta={np.degrees(theta[i]):.1f}°, phi={np.degrees(phi[i]):.2f}°, vel={vel_flat[i]} m/s")

    # Radar location
    lat_radar = radar.latitude['data'][0]
    lon_radar = radar.longitude['data'][0]

    # Convert Joplin to x/y relative to radar
    az, back_az, dist = geod.inv(lon_radar, lat_radar, lon_city, lat_city)
    az_rad = np.radians(az)
    x_city = dist/1e3 * np.sin(az_rad)
    y_city = dist/1e3 * np.cos(az_rad)

    # --- Plot ---
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure(figsize=(8, 7))
    ax = plt.subplot(111)

    display.plot_ppi(
        'velocity',
        sweep=sweep,
        ax=ax,
        cmap='RdBu_r',
        colorbar_label='Radial Velocity (m/s)',
        title=f'0.5 deg tilt Radial Velocity {timestamp}',
        vmax=40,
        vmin=-40,
        mask_outside=False
    )

    # Plot Joplin
    ax.plot(x_city, y_city, marker='*', color='green', markersize=10)
    ax.legend()

    # Range rings and crosshair
    display.plot_range_rings([50, 100, 150, 200], ax=ax, col='gray', ls='--')
    display.plot_cross_hair(5, ax=ax)
    plt.axis("equal")
    # Set plot limits ±50 km around city
    buffer = 40.0  # 50 km in meters
    plt.xlim([x_city - buffer, x_city + buffer])
    plt.ylim([y_city - buffer, y_city + buffer])

    # Save figure
    out_file = os.path.join(output_folder, f"radar_{timestamp}.png")
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)  # close figure to save memory

    print(f"Saved {out_file}")


