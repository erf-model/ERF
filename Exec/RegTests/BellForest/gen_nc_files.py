"""
Generate bell-shaped NetCDF input files for the BellForest test case.

Domain: 2000 m x 2000 m, center at (1000, 1000).

Files produced (all use x/y coordinate variables in metres):
  terrain_bell.nc      -- Gaussian terrain, peak 100 m, sigma 400 m
  forest_height_bell.nc -- Gaussian canopy height, peak 20 m, sigma 250 m
  forest_lai_bell.nc    -- Gaussian LAI, peak 2.0, sigma 250 m

Usage:
  python gen_nc_files.py
  (ncgen is called automatically; requires ncgen on PATH)
"""

import math
import os
import subprocess
import sys

NCGEN = "/opt/cray/pe/netcdf/4.9.0.3/bin/ncgen"

# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------
NX = 41          # number of x points  (matches 40-cell domain nodes: 0..2000 step 50)
NY = 41          # number of y points
DX = 50.0        # metres
DY = 50.0

cx = (NX - 1) * DX / 2.0   # 1000 m
cy = (NY - 1) * DY / 2.0   # 1000 m

x_vals = [i * DX for i in range(NX)]
y_vals = [j * DY for j in range(NY)]


def gaussian_2d(x, y, amp, sigma):
    """Gaussian bell centred on (cx, cy)."""
    r2 = (x - cx)**2 + (y - cy)**2
    return amp * math.exp(-r2 / (2.0 * sigma**2))


def make_field(amp, sigma):
    """Return [NY][NX] list (row-major, y outer) of Gaussian values."""
    rows = []
    for j in range(NY):
        for i in range(NX):
            rows.append(gaussian_2d(x_vals[i], y_vals[j], amp, sigma))
    return rows


def fmt_vals(vals, per_line=8):
    """Format a list of floats as a comma-separated CDL data string."""
    lines = []
    for start in range(0, len(vals), per_line):
        chunk = vals[start:start + per_line]
        lines.append(", ".join("{:.6f}".format(v) for v in chunk))
    return ",\n            ".join(lines)


def write_cdl(filename, varname, long_name, units, amp, sigma):
    """Write a CDL file for one 2-D Gaussian field."""
    data = make_field(amp, sigma)
    cdl = """\
netcdf {stem} {{
dimensions:
    y = {NY} ;
    x = {NX} ;
variables:
    double x(x) ;
        x:units = "m" ;
        x:long_name = "x coordinate" ;
    double y(y) ;
        y:units = "m" ;
        y:long_name = "y coordinate" ;
    double {var}(y, x) ;
        {var}:units = "{units}" ;
        {var}:long_name = "{long_name}" ;
        {var}:amplitude = {amp} ;
        {var}:sigma_m   = {sigma} ;
        {var}:cx_m      = {cx} ;
        {var}:cy_m      = {cy} ;
data:
    x = {x_vals} ;
    y = {y_vals} ;
    {var} =
            {data} ;
}}
""".format(
        stem=os.path.splitext(filename)[0],
        NY=NY, NX=NX,
        var=varname,
        long_name=long_name,
        units=units,
        amp=amp, sigma=sigma, cx=cx, cy=cy,
        x_vals=", ".join("{:.1f}".format(v) for v in x_vals),
        y_vals=", ".join("{:.1f}".format(v) for v in y_vals),
        data=fmt_vals(data),
    )
    cdl_file = filename.replace(".nc", ".cdl")
    with open(cdl_file, "w") as f:
        f.write(cdl)
    return cdl_file


def ncgen(cdl_file, nc_file):
    cmd = [NCGEN, "-o", nc_file, cdl_file]
    result = subprocess.call(cmd)
    if result != 0:
        sys.exit("ncgen failed for " + cdl_file)
    print("  wrote " + nc_file)


files = [
    # (nc filename,           CDL varname, long name,                units,  amp,   sigma)
    ("terrain_bell.nc",      "height",    "Gaussian bell terrain",   "m",   100.0, 400.0),
    ("forest_height_bell.nc","height",    "Gaussian canopy height",  "m",    20.0, 250.0),
    ("forest_lai_bell.nc",   "LAI",       "Gaussian leaf area index","m2/m2", 2.0, 250.0),
]

print("Generating NetCDF files in", os.getcwd())
for nc_file, varname, long_name, units, amp, sigma in files:
    cdl_file = write_cdl(nc_file, varname, long_name, units, amp, sigma)
    ncgen(cdl_file, nc_file)

print("Done.  CDL source files kept for inspection.")
