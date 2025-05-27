import pygrib
import numpy as np
import struct
from pyproj import Proj
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d
from IO import *
from Plot_1D import plot_1d
from Download_ERA5Data import *

const_g = 9.81

def p_sat(temp):
    tC = temp - 273.15  # Convert temperature from Kelvin to Celsius

    # Create masks for conditions
    mask_positive = tC > 0.0
    mask_negative = ~mask_positive

    # Initialize ps with zeros (same shape as temp)
    ps = np.zeros_like(temp)

    # Compute ps for tC > 0
    ps[mask_positive] = 6.112 * np.exp(17.62 * tC[mask_positive] / (tC[mask_positive] + 243.12))

    # Compute ps for tC <= 0
    ps[mask_negative] = 6.112 * np.exp(22.46 * tC[mask_negative] / (tC[mask_negative] + 272.62))

    return ps

def ReadERA5_3DData(file_path, is_IC):
    # Open the GRIB2 file
    pressure_levels = []
    ght_3d_hr3 = []
    temp_3d_hr3 = []
    uvel_3d_hr3 = []
    vvel_3d_hr3 = []
    wvel_3d_hr3 = []
    qv_3d_hr3 = []
    qc_3d_hr3 = []
    qr_3d_hr3 = []
    rh_3d_hr3 = []
    theta_3d_hr3 = []
    vort_3d_hr3 = []
    pressure_3d_hr3 = []
    lats, lons = None, None  # To store latitude and longitude grids

    printed_time = False

    datetime_str = ""

    with pygrib.open(file_path) as grbs:
        for grb in grbs:

            if not printed_time:
                year = grb.year
                month = grb.month
                day = grb.day
                hour = grb.hour
                minute = grb.minute if hasattr(grb, 'minute') else 0
                print(f"Date: {year}-{month:02d}-{day:02d}, Time: {hour:02d}:{minute:02d} UTC")
                datetime_str = f"{year:04d}_{month:02d}_{day:02d}_{hour:02d}_{minute:02d}"
                print(f"Datetime string: {datetime_str}")
                printed_time = True

            #print(f"Variable: {grb.name}, Level: {grb.level}, Units: {grb.parameterUnits}")
            if "Temperature" in grb.name:
                # Append temperature values
                temp_3d_hr3.append(grb.values)

                # Append pressure level
                pressure_levels.append(grb.level)

            if "Geopotential" in grb.name:
                ght_3d_hr3.append(grb.values)

            if "Potential temperature" in grb.name:
                theta_3d_hr3.append(grb.values)

            if "Pressure" in grb.name:
                pressure_3d_hr3.append(grb.values)

            if "U component of wind" in grb.name:
                uvel_3d_hr3.append(grb.values)

            if "V component of wind" in grb.name:
                vvel_3d_hr3.append(grb.values)

            if "Vertical velocity" in grb.name:
                wvel_3d_hr3.append(grb.values)

            if "Specific humidity" in grb.name:
                qv_3d_hr3.append(grb.values)

            if "Specific cloud liquid water content" in grb.name:
                qc_3d_hr3.append(grb.values)

            if "Specific rain water content" in grb.name:
                qr_3d_hr3.append(grb.values)

            if "Relative humidity" in grb.name:
                rh_3d_hr3.append(grb.values)

            if "Vorticity (relative)" in grb.name:
                vort_3d_hr3.append(grb.values)

            # Retrieve latitude and longitude grids (once)
            if lats is None or lons is None:
                lats, lons = grb.latlons()

    # Stack into a 3D array (level, lat, lon)
    ght_3d_hr3 = np.stack(ght_3d_hr3, axis=0)
    uvel_3d_hr3 = np.stack(uvel_3d_hr3, axis=0)
    vvel_3d_hr3 = np.stack(vvel_3d_hr3, axis=0)
    wvel_3d_hr3 = np.stack(wvel_3d_hr3, axis=0)
    #theta_3d_hr3 = np.stack(theta_3d_hr3, axis=0)
    qv_3d_hr3 = np.stack(qv_3d_hr3, axis=0)
    qc_3d_hr3 = np.stack(qc_3d_hr3, axis=0)
    qr_3d_hr3 = np.stack(qr_3d_hr3, axis=0)
    rh_3d_hr3 = np.stack(rh_3d_hr3, axis=0)
    temp_3d_hr3 = np.stack(temp_3d_hr3, axis=0)
    vort_3d_hr3 = np.stack(vort_3d_hr3, axis=0)

    #pressure_3d_hr3 = np.stack(pressure_3d_hr3, axis=0)
    # Get the size of each dimension
    dim1, dim2, dim3 = ght_3d_hr3.shape

    # Print the sizes
    print(f"Size of dimension 1: {dim1}")
    print(f"Size of dimension 2: {dim2}")
    print(f"Size of dimension 3: {dim3}")

    # Convert pressure levels to numpy array for indexing
    pressure_levels = np.array(pressure_levels)
    nz = 37

    print("The number of lats and lons are levels are %d, %d, %d"%(lats.shape[0], lats.shape[1], nz));


    # Extract unique latitude and longitude values
    unique_lats = np.unique(lats[:, 0])  # Take the first column for unique latitudes
    unique_lons = np.unique(lons[0, :])  # Take the first row for unique longitudes

    print("Min max lat lons are ", unique_lats[0], unique_lats[-1], unique_lons[0], unique_lons[-1]);

    nlats = len(unique_lats)
    nlons = len(unique_lons)

    lat_start = 0
    lat_end = nlats
    lon_start = 0
    lon_end = nlons

    domain_lats = unique_lats[:]
    domain_lons = unique_lons[:]

    print("The min max are",(lat_start, lat_end, lon_start, lon_end));

    nx = domain_lats.shape[0]
    ny = domain_lons.shape[0]

    print("nx and ny here are ", nx, ny)

    z_grid = np.zeros((nx, ny, nz))
    rhod_3d = np.zeros((nx, ny, nz))
    uvel_3d = np.zeros((nx, ny, nz))
    vvel_3d = np.zeros((nx, ny, nz))
    wvel_3d = np.zeros((nx, ny, nz))
    theta_3d = np.zeros((nx, ny, nz))
    qv_3d = np.zeros((nx, ny, nz))
    qc_3d = np.zeros((nx, ny, nz))
    qr_3d = np.zeros((nx, ny, nz))
    rh_3d = np.zeros((nx, ny, nz))
    temp_3d = np.zeros((nx, ny, nz))
    qsat_3d = np.zeros((nx, ny, nz))

    velocity = np.zeros((nx, ny, nz, 3))

    vort_3d = np.zeros((nx, ny, nz))
    pressure_3d = np.zeros((nx, ny, nz))
    theta_3d = np.zeros((nx, ny, nz))


    # Create meshgrid
    x_grid, y_grid = np.meshgrid(domain_lons, domain_lats)
    lon_grid, lat_grid = np.meshgrid(domain_lons, domain_lats)

    utm_zone = calculate_utm_zone(-90.50399843097345)
    utm_proj = Proj(proj='utm', zone=utm_zone, ellps='WGS84')

    # Convert the entire grid to UTM
    x_grid, y_grid = utm_proj(lon_grid, lat_grid)

    k_to_delete = []

    print("size is ", len(qv_3d_hr3))

    pressure_typical = np.loadtxt("./TypicalAtmosphereData/pressure_vs_z_actual.txt")
    pressure_interp_func = interp1d(pressure_typical[:,1], pressure_typical[:,0], kind='linear', fill_value="extrapolate")

    # Find the index of the desired pressure level
    for k in np.arange(36, -1, -1):

        # Extract temperature at the desired pressure level
        ght_at_lev = ght_3d_hr3[k]
        #pressure_at_lev = pressure_3d_hr3[k]
        uvel_at_lev = uvel_3d_hr3[k]
        vvel_at_lev = vvel_3d_hr3[k]
        if(k >= len(wvel_3d_hr3)):
            wvel_at_lev = 0.0
        else:
            wvel_at_lev = wvel_3d_hr3[k]
        #theta_at_lev = theta_3d_hr3[k]
        rh_at_lev = rh_3d_hr3[k]
        if(k >= len(qv_3d_hr3)):
            qv_at_lev = 0.0
        else:
            qv_at_lev = qv_3d_hr3[k]
        if(k >= len(qc_3d_hr3)):
            qc_at_lev = 0.0
        else:
            qc_at_lev = qc_3d_hr3[k]
        if(k >= len(qr_3d_hr3)):
            qr_at_lev = 0.0
        else:
            qr_at_lev = qr_3d_hr3[k]
        temp_at_lev = temp_3d_hr3[k]
        #if(k >= len(theta_3d_hr3)):
        #    theta_at_lev = 0.0
        #else:
        #    theta_at_lev = theta_3d_hr3[k]
        if(k >= len(vort_3d_hr3)):
            vort_at_lev = 0.0
        else:
            vort_at_lev = vort_3d_hr3[k]

        # Print temperature for the specified level using nested loops
        #for i in range(lats.shape[0]):  # Latitude index
            #for j in range(lats.shape[1]):  # Longitude index
        #lat = lats[i, j]
        #lon = lons[i, j]
        temp_3d[:, :, k] = temp_at_lev
        z_grid[:,:,k] = ght_at_lev/const_g
        #pressure_3d[:,:,k] = pressure_at_lev

        uvel_3d[:, :, k] = uvel_at_lev
        vvel_3d[:, :, k] = vvel_at_lev
        wvel_3d[:, :, k] = wvel_at_lev
        #theta_3d[:, :, k] = theta_at_lev
        rh_3d[:, :, k] = rh_at_lev
        qv_3d[:, :, k] = qv_at_lev
        qc_3d[:, :, k] = qc_at_lev
        qr_3d[:, :, k] = qr_at_lev
        vort_3d[:, :, k] = vort_at_lev

        qv_3d[:, :, k] = qv_at_lev
        rh_val = rh_at_lev

        ps = p_sat(temp_3d[:, :, k]);
        pv = rh_val/100.0*ps;

        print("Max val is ", np.max(z_grid[:,:,k]),  )



        #pressure_3d[:, :, k] = (temp_3d[:, :, k]/theta_3d[:, :, k])**(1004.5/287.0)*1000.0
        #pressure_3d[:, :, k] = 0.622*pv/qv_3d[:, :, k] + pv
        #pressure_3d[:, :, k] = 1000.0*np.exp(-const_g*(np.mean(z_grid[:,:,k])-0.0)/(287*temp_3d[:, :, k]*(1.0+1.6*qv_3d[:, :, k])))

        # Assuming quantities at surface is same as the first cell
        if(k==36):
            pressure_3d[:, :, k] = 1000.0 - 1000.0/(287*temp_3d[:, :, k]*(1.0+1.6*qv_3d[:, :, k]))*const_g*z_grid[:,:,k]
        else:
            pressure_3d[:, :, k] = pressure_3d[:, :, k+1] - pressure_3d[:, :, k+1]/(287*temp_3d[:, :, k+1]*(1.0+1.6*qv_3d[:, :, k+1]))*const_g*(z_grid[:,:,k]-z_grid[:,:,k+1])

        qsat_3d[:,:,k] = 0.622*ps/(pressure_3d[:, :, k]-ps)

        #pressure_typical_here = pressure_interp_func(np.mean(z_grid[:,:,k]));
        #indices = np.argwhere( (pressure_3d[:,:,k] <= 0.9*pressure_typical_here) | (pressure_3d[:,:,k] >= 1.02*pressure_typical_here))

        #for index in indices:
            #row, col = index
            #pressure_3d[row, col,k] = pressure_typical_here

        rhod_3d[:,:,k] = pressure_3d[:, :, k]*100.0/(287.0*temp_3d[:, :, k]*(1.0 + 1.6*qv_3d[:, :, k]))

        #condition = np.any((pressure_3d[:,:,k] <= 1e-3))
        #if(condition):
            #print("k val is ", k, np.min(pressure_3d[:,:,k]), np.max(pressure_3d[:,:,k]), np.min(temp_3d[:, :, k]), np.min(rh_val), np.max(rh_val), np.max(qv_3d[:, :, k]))
            #k_to_delete.append(k)

        theta_3d[:,:,k] = temp_3d[:, :, k]*(1000.0/pressure_3d[:, :, k])**(287.0/1004.5)


        # Find indices of elements that are zero or less
        #indices = np.argwhere(qv_3d[:, :, k] <= 0)
        indices = np.argwhere(rh_val <= 0)

        # Print indices and values
        for index in indices:
            row, col = index
            value = rh_val[row, col]
            #print(f"Element at index ({row}, {col}, {k}) is zero or less: {value}")

        velocity[:,:,k,0] = uvel_at_lev
        velocity[:,:,k,1] = vvel_at_lev
        velocity[:,:,k,2] = 0.0


        #print(f"Lat and lon are: {lat_grid[0,0]:.2f}, {lon_grid[0,0]:.2f}")
        #print(f"Temperature: {temp_3d[0,0,k]:.2f} K, Pressure: {pressure_3d[0,0,k]:.2f}, Geo height : {z_grid[0,0,k]:.2f} ")


    scalars = {
         "density": rhod_3d,
         "uvel": uvel_3d,
         "vvel": vvel_3d,
         "wvel": wvel_3d,
         "theta": theta_3d,
         "qv": qv_3d,
         "qc": qc_3d,
         "qr": qr_3d,
         "rh": rh_3d,
         "temperature": temp_3d,
         "vorticity": vort_3d,
         "pressure": pressure_3d,
         "qsat": qsat_3d,
    }


    dir_path = "Images"
    os.makedirs(dir_path, exist_ok=True)

    dir_path = "Output"
    os.makedirs(dir_path, exist_ok=True)

    output_vtk = "./Output/ERA5_" + datetime_str + ".vtk"

    output_binary = "./Output/ERF_IC_" + datetime_str + ".bin"

    write_binary_vtk_structured_grid(output_vtk, x_grid, y_grid, z_grid,
                                     nz, k_to_delete, True,
                                     scalars, velocity)

    write_binary_vtk_cartesian(output_binary, domain_lats, domain_lons,
                               x_grid, y_grid, z_grid,
                               nx, ny, nz, k_to_delete, scalars)

    scalars_for_ERF = {
         "theta": theta_3d,
         "uvel": uvel_3d,
         "vvel": vvel_3d,
    }

    plot_1d(temp_3d, pressure_3d, theta_3d, qv_3d, qsat_3d, z_grid, k_to_delete)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 download_era5.py <input_filename>")
        sys.exit(1)

input_filename = sys.argv[1]
filename  = download_era5_data(input_filename);
print("Filename is ", filename);

file_paths = [filename]

# Loop through each file and process it
count = 0;
is_IC = True
for file_path in file_paths:
    print(f"Processing file: {file_path}")
    #ReadERA5_SurfaceData("Laura_SurfaceData.grib", is_IC)
    ReadERA5_3DData(file_path, is_IC)
    is_IC = False
