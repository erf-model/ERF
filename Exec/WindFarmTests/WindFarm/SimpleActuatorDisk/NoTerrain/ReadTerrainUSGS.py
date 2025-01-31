import math
import numpy as np
import sys
import os
import rasterio
from rasterio.transform import xy
from rasterio.warp import transform_bounds
from rasterio.transform import rowcol

# Check if two filename arguments are provided
if len(sys.argv) != 3:
    print("Usage: python3 <script_name> <filename1> <filename2>")
    sys.exit(1)

# Get the filenames from the command-line arguments
filename1 = sys.argv[1]
filename2 = sys.argv[2]

# Check if both files exist, throw an error if not
if not os.path.isfile(filename1):
    raise FileNotFoundError(f"Error: The file '{filename1}' does not exist.")
if not os.path.isfile(filename2):
    raise FileNotFoundError(f"Error: The file '{filename2}' does not exist.")


# Path to the GeoTIFF file
geotiff_file = filename1

domain_bounds = np.loadtxt(filename2)

# Define domain bounds
domain_lon_min = domain_bounds[0]
domain_lon_max = domain_bounds[1]
domain_lat_min = domain_bounds[2]
domain_lat_max = domain_bounds[3]

# Constants
rad_earth = 6371000.0  # Radius of Earth in meters
M_PI = math.pi
dist_per_deg_lat = 6371000.0*2.0*M_PI/(2.0*180.0)

row = -1
col = -1

domain_row_min = 0
domain_row_max = 0
domain_col_min = 0
domain_col_max = 0

# Open the GeoTIFF file
with rasterio.open(geotiff_file) as src:
    # Get the affine transform of the raster
    transform = src.transform

    # Calculate the row corresponding to the given latitude
    domain_row_min = int((src.bounds.top - domain_lat_max) / abs(transform[4]))  # transform[4] is the pixel height (negative)
    domain_row_max = int((src.bounds.top - domain_lat_min) / abs(transform[4]))  # transform[4] is the pixel height (negative)
    domain_col_min = int((domain_lon_min - transform[2]) / transform[0])  # transform[2] is the x-coordinate of the top-left corner
    domain_col_max = int((domain_lon_max - transform[2]) / transform[0])  # transform[2] is the x-coordinate of the top-left corner

    print("src.bounds_top are %0.15g %0.15g\n"%(src.bounds.top, transform[4]));
    print("The row and cols are %d, %d, %d, %d"%(domain_row_min, domain_row_max, domain_col_min, domain_col_max));

     # Convert row and column back to latitude and longitude
    lon_min = transform[2] + domain_col_min * transform[0]
    lon_max = transform[2] + domain_col_max * transform[0]

    lat_min = src.bounds.top + domain_row_max * transform[4]
    lat_max = src.bounds.top + domain_row_min * transform[4]

    print("The lon min and max are %0.15g, %0.15g"%(lon_min, lon_max));
    print("The lat min and max are %0.15g, %0.15g"%(lat_min, lat_max));


    print("Values are %.15g, %0.15g, %0.15g, %0.15g\n"%(src.bounds.top,domain_row_max,transform[4], lat_min));
# Open the GeoTIFF file
with rasterio.open(geotiff_file) as src:

    # Read the raster shape (rows, columns)
    rows, cols = src.shape
    transform = src.transform

    # Print CRS and file bounds for context
    print(f"CRS: {src.crs}")
    print(f"Bounds: {src.bounds}")
    print(f"Raster Shape (rows, cols): {rows}, {cols}\n")

    if(domain_bounds[0] <= src.bounds.left or domain_bounds[0] >= src.bounds.right):
        print("The longitude min (the first entry) specified in %s is %0.15g,"
               " which is not within the data longitude bounds (%0.15g, %0.15g). Exiting....\n"%
               (filename2, domain_bounds[0], src.bounds.left,src.bounds.right));
        sys.exit();

    if(domain_bounds[1] <= src.bounds.left or domain_bounds[1] >= src.bounds.right):
        print("The longitude max (the second entry) specified in %s is %0.15g,"
               " which is not within the data longitude bounds (%0.15g, %0.15g). Exiting....\n"%
               (filename2, domain_bounds[1], src.bounds.left,src.bounds.right));
        sys.exit();

    if(domain_bounds[2] <= src.bounds.bottom or domain_bounds[2] >= src.bounds.top):
        print("The latitude min (the third entry) specified in %s is %0.15g,"
               " which is not within the data longitude bounds (%0.15g, %0.15g). Exiting....\n"%
               (filename2, domain_bounds[2], src.bounds.bottom,src.bounds.top));
        sys.exit();

    if(domain_bounds[3] <= src.bounds.bottom or domain_bounds[3] >= src.bounds.top):
        print("The latitude max (the fourth entry) specified in %s is %0.15g,"
               " which is not within the data longitude bounds (%0.15g, %0.15g). Exiting....\n"%
               (filename2, domain_bounds[3], src.bounds.bottom,src.bounds.top));
        sys.exit();

    # Get the bounds in the source CRS
    bounds = src.bounds

    # Transform bounds to latitude and longitude (EPSG:4326)
    lat_lon_bounds = transform_bounds(src.crs, "EPSG:4326", bounds.left, bounds.bottom, bounds.right, bounds.top)

    print("Data lat min %0.15g\n"%bounds.bottom);
    print("Data lat max %0.15g\n"%bounds.top);
    print("Data lon min %0.15g\n"%bounds.left);
    print("Data lon max %0.15g\n"%bounds.right);

    lat_min = lat_min*M_PI / 180.0
    lon_min = lon_min*M_PI / 180.0

    # Read the first band of elevation data
    elevation_data = src.read(1)

    xloc = 0;
    yloc = []  # To store y-coordinates

    vtk_file = open("terrain_mesh.vtk",'w')

    vtk_file.write("# vtk DataFile Version 2.0\n")
    vtk_file.write("Structured Grid Example\n")
    vtk_file.write("ASCII\n\n")

    # Write the structured grid header
    vtk_file.write(f"DATASET STRUCTURED_GRID\n")
    nskip = 1
    rows_in_output = int((domain_row_max-domain_row_min)/nskip);
    cols_in_output = int((domain_col_max-domain_col_min)/nskip);
    vtk_file.write(f"DIMENSIONS {cols_in_output} {rows_in_output} 1\n")  # Adjust for structured grid dimensions
    vtk_file.write(f"POINTS {rows_in_output*cols_in_output} float\n")  # Adjust for structured grid dimensions

    file_for_erf = open("ERF_terrain_file.txt",'w')
    file_for_erf.write("%0.15g %0.15g\n"%(lon_min*180/M_PI,lat_min*180/M_PI))
    file_for_erf.write("%d %d\n"%(domain_col_max-domain_col_min,domain_row_max-domain_row_min));

    for row in range(domain_row_max, domain_row_min, -1):  # Sample every 10% of rows
        print("Doing %d of %d"%(row-domain_row_min, domain_row_max-domain_row_min))
        for col in range(domain_col_min, domain_col_max, nskip):  # Sample every 10% of columns

            # Read the value at this pixel location
            value = elevation_data[row, col]

            # Convert row and column back to latitude and longitude
            lon = transform[2] + col * transform[0]  # lon at the given column
            lat = src.bounds.top + row * transform[4]  # lat at the given row

            # Get the elevation value
            elevation = elevation_data[row, col]

            if(elevation==0.0):
                print("Elevation is zero...cannot be\n");
                sys.exit()

            # Conversion and calculations
            lat = lat * M_PI / 180.0
            lon = lon * M_PI / 180.0
            delta_lat = lat - lat_min
            delta_lon = lon - lon_min

            term1 = math.pow(math.sin(delta_lat / 2.0), 2)
            term2 = math.cos(lat) * math.cos(lat_min) * math.pow(math.sin(delta_lon / 2.0), 2)
            dist = 2.0 * rad_earth * math.asin(math.sqrt(term1 + term2))
            dy = (lat - lat_min) * dist_per_deg_lat * 180.0 / M_PI
            fac = math.pow(dist, 2) - math.pow(dy, 2)
            if(dist < dy):
                if(abs(dist-dy) > 1e-10):
                    print("Error in calculation dist cannot be less than dy %0.15g, %0.15g"%(dist, dy))
                    sys.exit()

            if(abs(dist-dy) < 1e-10):
                fac = 0.0

            dx = math.sqrt(fac)

            xloc = dx
            file_for_erf.write("%0.15g %0.15g %0.15g\n"%(xloc,dy,elevation));

            vtk_file.write(f"{xloc} {dy} {elevation}\n")
            if(col==domain_col_min):
                yloc.append(dy)

    file_for_erf.close()

# Write the elevation as a scalar field
    vtk_file.write(f"\nPOINT_DATA {rows_in_output * cols_in_output}\n")
    vtk_file.write("SCALARS Elevation float 1\n")
    vtk_file.write("LOOKUP_TABLE default\n")
    for row in range(domain_row_max, domain_row_min, -1):
        for col in range(domain_col_min, domain_col_max, nskip):
            elevation_value = elevation_data[row, col]  # Elevation for the current point
            vtk_file.write(f"{elevation_value}\n")

file_for_erf.close()

