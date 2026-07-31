import sys
import os
import argparse
from pyproj import Transformer
import numpy as np
import pyvista as pv
import cartopy.feature as cfeature

def read_user_input(filename):
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            if ':' in line:
                key, value = line.strip().split(':', 1)
                key = key.strip().lower()
                value = value.strip()
                if key == 'area':
                    data[key] = [float(x) for x in value.split(',')]
                elif key == 'time':
                    data[key] = value
                else:
                    data[key] = int(value)
    return data

def CreateLCCMapping(area):
    lat1 = area[2]
    lat2 = area[0]
    lon1 = area[1]
    lon2 = area[3]

    delta = lat2 - lat1
    lon0 = (lon1 + lon2) / 2
    lat0 = (lat1 + lat2) / 2

    lat_1 = lat1 + delta/6
    lat_2 = lat2 - delta/6

    lambert_conformal = (
        f"+proj=lcc +lat_1={lat_1:.6f} +lat_2={lat_2:.6f} "
        f"+lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs"
    )
    return lambert_conformal

def WriteUSMapLegacyVTK(filename, all_points, lines_flat):
    """
    Write legacy ASCII VTK PolyData file (works with VisIt and ParaView).
    """
    n_points = all_points.shape[0]
    
    # Count number of lines
    n_lines = 0
    idx = 0
    while idx < len(lines_flat):
        n = lines_flat[idx]
        n_lines += 1
        idx += n + 1

    with open(filename, "w") as f:
        # Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("US state borders in Lambert projection\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        # Points
        f.write(f"POINTS {n_points} float\n")
        for x, y, z in all_points:
            f.write(f"{x} {y} {z}\n")

        # Lines
        f.write(f"LINES {n_lines} {len(lines_flat)}\n")
        idx = 0
        while idx < len(lines_flat):
            n = lines_flat[idx]
            indices = lines_flat[idx+1:idx+1+n]
            indices_str = " ".join(str(i) for i in indices)
            f.write(f"{n} {indices_str}\n")
            idx += n + 1
    print(f"Saved legacy VTK: {filename}")



def WriteUSMapVTKFile(area, filename="USMap_LambertProj.vtk"):
    # Define bounds (clip to continental US or user area if you prefer)
    lon_min, lon_max = -125, -50
    lat_min, lat_max = 24, 50

    # Get US state borders from cartopy
    states_feature = cfeature.STATES.with_scale('50m')
    geoms = list(states_feature.geometries())

    # Create transformer FROM geographic (lon/lat) TO Lambert
    lambert_conformal = CreateLCCMapping(area)
    transformer = Transformer.from_crs("EPSG:4326", lambert_conformal, always_xy=True)

    all_points = []
    lines_flat = []
    pt_id = 0

    for geom in geoms:
        try:
            # MultiLineString
            for line in geom.boundary.geoms:
                coords = [(x, y) for x, y in line.coords if lon_min <= x <= lon_max and lat_min <= y <= lat_max]
                if len(coords) < 2:
                    continue
                # Project coordinates to Lambert
                x_proj, y_proj = transformer.transform([c[0] for c in coords],
                                                       [c[1] for c in coords])
                coords_proj = list(zip(x_proj, y_proj))
                n = len(coords_proj)
                all_points.extend([[x, y, 0] for x, y in coords_proj])
                lines_flat.extend([n] + list(range(pt_id, pt_id+n)))
                pt_id += n
        except AttributeError:
            # Single LineString
            coords = [(x, y) for x, y in geom.boundary.coords if lon_min <= x <= lon_max and lat_min <= y <= lat_max]
            if len(coords) < 2:
                continue
            # Project
            x_proj, y_proj = transformer.transform([c[0] for c in coords],
                                                   [c[1] for c in coords])
            coords_proj = list(zip(x_proj, y_proj))
            n = len(coords_proj)
            all_points.extend([[x, y, 0] for x, y in coords_proj])
            lines_flat.extend([n] + list(range(pt_id, pt_id+n)))
            pt_id += n

    # Convert to numpy arrays
    all_points = np.array(all_points)
    lines_flat = np.array(lines_flat, dtype=np.int64)

    # Create PolyData
    polydata = pv.PolyData()
    polydata.points = all_points
    polydata.lines = lines_flat

    # Save to VTK
    polydata.save(filename)
    WriteUSMapLegacyVTK("USMap_LambertProj.vtk", all_points, lines_flat)

    print(f"Saved {filename}")
    return lambert_conformal


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write USMap in Lambert projection coordinates to ASCII VTK")
    parser.add_argument("input_filename", help="Input file with area/time/etc.")
    args = parser.parse_args()

    user_inputs = read_user_input(args.input_filename)
    print("User inputs:", user_inputs)
    area = user_inputs.get("area", None)

    lambert_conformal = WriteUSMapVTKFile(area)

