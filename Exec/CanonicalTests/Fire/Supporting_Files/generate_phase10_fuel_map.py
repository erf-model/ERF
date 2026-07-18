#!/usr/bin/env python3
"""
Generate spatial fuel map files for Phase 10 regression tests.

Produces ESRI ASCII grid format fuel maps for use by the ERF fire simulation.
Files are compatible with the read_ascii_fuel_map() function in ERF_FuelMap.H.

ESRI ASCII Grid Format:
  - Header: ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value
  - Data: rows ordered from north (highest y) to south (lowest y)
  - Fire domain then reverses row order to match domain coordinates (south to north)

Domain Parameters:
  - Grid ratio: 5 (fire cells are 5x finer than atmospheric cells)
  - Atmospheric cells: 40 x 40 at 50 m resolution
  - Fire cells: 200 x 200 at 10 m resolution
  - Total domain: 2000 m x 2000 m

Generated Files:
  1. fuel_map_phase10_200x200.asc
     - Two-zone map with fuel model 1 (FM1, short grass) in southern half (y < 1000 m)
     - Fuel model 4 (FM4, chaparral) in northern half (y >= 1000 m)
     - Used by spatial_fuel and blending regression tests

  2. fuel_map_phase10_firebreak.asc
     - Uniform fuel model 1 (FM1) across entire domain
     - Used by firebreak regression test
"""

import sys
import os


def create_esri_ascii_header(ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value):
    """
    Create ESRI ASCII grid header lines.
    
    Args:
        ncols: Number of columns (cells in x direction)
        nrows: Number of rows (cells in y direction)
        xllcorner: X coordinate of lower-left corner [m]
        yllcorner: Y coordinate of lower-left corner [m]
        cellsize: Cell size [m]
        nodata_value: Value representing missing/invalid data
    
    Returns:
        List of header lines (without newlines)
    """
    header = [
        f"ncols {ncols}",
        f"nrows {nrows}",
        f"xllcorner {xllcorner}",
        f"yllcorner {yllcorner}",
        f"cellsize {cellsize}",
        f"nodata_value {nodata_value}"
    ]
    return header


def generate_fuel_map_two_zone(ncols=200, nrows=200, boundary_y=1000.0, 
                                cellsize=10.0, fm_south=1, fm_north=4):
    """
    Generate a two-zone fuel map: FM1 (south) and FM4 (north) separated at y=boundary_y.
    
    In ESRI format, rows are stored from north (highest y) to south (lowest y).
    Row 0 in file corresponds to domain y = 2000 m (northernmost).
    Row nrows-1 in file corresponds to domain y = 0 m (southernmost).
    
    The fire solver will reverse row order internally to match domain coordinates.
    
    Args:
        ncols: Number of columns (fire grid in x) [default 200]
        nrows: Number of rows (fire grid in y) [default 200]
        boundary_y: Y-coordinate boundary between FM zones [m, default 1000]
        cellsize: Fire grid cell size [m, default 10.0]
        fm_south: Fuel model code for southern zone [default 1]
        fm_north: Fuel model code for northern zone [default 4]
    
    Returns:
        Tuple: (header_lines, data_rows) where data_rows is list of lists
    """
    # Domain extent
    prob_lo_y = 0.0
    prob_hi_y = nrows * cellsize
    
    header = create_esri_ascii_header(ncols, nrows, 0.0, prob_lo_y, cellsize, -9999)
    
    data_rows = []
    
    # Iterate from north (highest y) to south (lowest y) as ESRI format requires
    for file_row in range(nrows):
        # Map file row to domain y coordinate
        # file_row = 0 => y = prob_hi_y (northernmost)
        # file_row = nrows-1 => y = prob_lo_y (southernmost)
        domain_y = prob_hi_y - (file_row + 0.5) * cellsize
        
        # Determine fuel code for this row based on boundary
        fuel_code = fm_north if domain_y >= boundary_y else fm_south
        
        # Generate row of fuel codes (all same code for this row)
        row = [fuel_code] * ncols
        data_rows.append(row)
    
    return header, data_rows


def generate_fuel_map_uniform(ncols=200, nrows=200, cellsize=10.0, fm_code=1):
    """
    Generate a uniform fuel map with single fuel model everywhere.
    
    Args:
        ncols: Number of columns (fire grid in x) [default 200]
        nrows: Number of rows (fire grid in y) [default 200]
        cellsize: Fire grid cell size [m, default 10.0]
        fm_code: Fuel model code to use throughout [default 1]
    
    Returns:
        Tuple: (header_lines, data_rows) where data_rows is list of lists
    """
    prob_lo_y = 0.0
    
    header = create_esri_ascii_header(ncols, nrows, 0.0, prob_lo_y, cellsize, -9999)
    
    # All rows identical: uniform fuel code
    data_rows = [[fm_code] * ncols for _ in range(nrows)]
    
    return header, data_rows


def write_esri_ascii_file(filename, header, data_rows):
    """
    Write fuel map to ESRI ASCII grid file.
    
    Args:
        filename: Path to output file
        header: List of header lines (strings)
        data_rows: List of lists, each containing fuel codes for a row
    """
    with open(filename, 'w') as f:
        # Write header
        for line in header:
            f.write(line + "\n")
        
        # Write data rows
        for row in data_rows:
            f.write(" ".join(str(code) for code in row) + "\n")


def main():
    """Generate fuel map files for Phase 10 tests."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # File 1: Two-zone map (FM1 south, FM4 north)
    header1, data1 = generate_fuel_map_two_zone(ncols=200, nrows=200, 
                                                boundary_y=1000.0, cellsize=10.0, 
                                                fm_south=1, fm_north=4)
    file1 = os.path.join(script_dir, "fuel_map_phase10_200x200.asc")
    write_esri_ascii_file(file1, header1, data1)
    print(f"Created: {file1}")
    
    # File 2: Uniform FM1 map
    header2, data2 = generate_fuel_map_uniform(ncols=200, nrows=200, 
                                               cellsize=10.0, fm_code=1)
    file2 = os.path.join(script_dir, "fuel_map_phase10_firebreak.asc")
    write_esri_ascii_file(file2, header2, data2)
    print(f"Created: {file2}")


if __name__ == "__main__":
    main()
