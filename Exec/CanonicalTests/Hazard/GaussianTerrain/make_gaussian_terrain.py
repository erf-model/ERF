#!/usr/bin/env python3
"""
Generate Gaussian terrain files in ERF's non-USGS terrain format.

Usage:
    python make_gaussian_terrain.py \
        --nx 32 --ny 32 \
        --xlo 0.0 --ylo 0.0 \
        --dx 125.0 --dy 125.0 \
        --height 200.0 \
        --sigma 600.0 \
        --output gaussian_hill.txt

The terrain file format is:
    nx
    ny
    x[0]
    ...
    x[nx-1]
    y[0]
    ...
    y[ny-1]
    z[0,0]
    z[0,1]
    ...
    z[nx-1,ny-1]

with one value per line and z stored contiguous in y
(outer loop over x/i, inner loop over y/j).

The Gaussian formula is:
    z(x,y) = H * exp(-((x-cx)^2 + (y-cy)^2) / (2*sigma^2))

where cx = xlo + nx*dx/2, cy = ylo + ny*dy/2 (domain centre).
"""

import argparse
import math


def generate_gaussian_terrain(nx, ny, xlo, ylo, dx, dy, height, sigma, output_file):
    """
    Generate a Gaussian terrain file in ERF's non-USGS terrain format.

    Parameters
    ----------
    nx, ny : int
        Number of grid points in x and y
    xlo, ylo : float
        Lower-left corner of domain
    dx, dy : float
        Grid spacing in x and y
    height : float
        Gaussian amplitude (positive for hill, negative for pit)
    sigma : float
        Gaussian width parameter (standard deviation)
    output_file : str
        Path to output file
    """
    # Domain center
    cx = xlo + 0.5 * nx * dx
    cy = ylo + 0.5 * ny * dy

    x_coords = [xlo + i * dx for i in range(nx)]
    y_coords = [ylo + j * dy for j in range(ny)]

    # Precompute 2*sigma^2
    two_sigma2 = 2.0 * sigma * sigma

    z_values = []
    for x in x_coords:
        for y in y_coords:
            r_squared = (x - cx) ** 2 + (y - cy) ** 2
            z_values.append(height * math.exp(-r_squared / two_sigma2))

    # Write one value per line: nx, ny, x-array, y-array, then z contiguous in y
    with open(output_file, 'w') as f:
        f.write(f"{nx}\n")
        f.write(f"{ny}\n")

        for x in x_coords:
            f.write(f"{x:.1f}\n")

        for y in y_coords:
            f.write(f"{y:.1f}\n")

        for z in z_values:
            f.write(f"{z:14.7f}\n")

    print(f"Written {output_file}")

    z_min = min(z_values)
    z_max = max(z_values)
    print(f"  Terrain stats: min={z_min:.3f}, max={z_max:.3f}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Gaussian terrain file in ERF non-USGS terrain format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--nx", type=int, default=32, help="Grid points in x")
    parser.add_argument("--ny", type=int, default=32, help="Grid points in y")
    parser.add_argument("--xlo", type=float, default=0.0, help="Lower x coordinate")
    parser.add_argument("--ylo", type=float, default=0.0, help="Lower y coordinate")
    parser.add_argument("--dx", type=float, default=125.0, help="Grid spacing in x")
    parser.add_argument("--dy", type=float, default=125.0, help="Grid spacing in y")
    parser.add_argument("--height", type=float, default=200.0, help="Gaussian height (+ = hill, - = pit)")
    parser.add_argument("--sigma", type=float, default=600.0, help="Gaussian width parameter")
    parser.add_argument("--output", type=str, required=True, help="Output file path")

    args = parser.parse_args()

    generate_gaussian_terrain(
        nx=args.nx,
        ny=args.ny,
        xlo=args.xlo,
        ylo=args.ylo,
        dx=args.dx,
        dy=args.dy,
        height=args.height,
        sigma=args.sigma,
        output_file=args.output
    )


if __name__ == "__main__":
    main()
