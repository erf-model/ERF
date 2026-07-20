#!/usr/bin/env python3
"""
Generate Gaussian terrain files in Askervein format.

Usage:
    python make_gaussian_terrain.py \
        --nx 32 --ny 32 \
        --xlo 0.0 --ylo 0.0 \
        --dx 125.0 --dy 125.0 \
        --height 200.0 \
        --sigma 600.0 \
        --output gaussian_hill.txt

The Gaussian formula is:
    z(x,y) = H * exp(-((x-cx)^2 + (y-cy)^2) / (2*sigma^2))

where cx = xlo + nx*dx/2, cy = ylo + ny*dy/2 (domain centre).
"""

import argparse
import math


def generate_gaussian_terrain(nx, ny, xlo, ylo, dx, dy, height, sigma, output_file):
    """
    Generate a Gaussian terrain file in Askervein format.

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

    # Precompute 2*sigma^2
    two_sigma2 = 2.0 * sigma * sigma

    # Write to file
    with open(output_file, 'w') as f:
        # Header: nx, ny
        f.write(f"{nx}. {ny}\n")
        # Domain bounds and spacing
        f.write(f"{xlo}. {ylo}\n")
        f.write(f"{dx}. {dy}\n")

        # Generate and write z values row by row
        for j in range(ny):
            for i in range(nx):
                # Node coordinates
                x = xlo + i * dx
                y = ylo + j * dy

                # Gaussian formula
                r_squared = (x - cx) ** 2 + (y - cy) ** 2
                z = height * math.exp(-r_squared / two_sigma2)

                f.write(f"{z:14.7f}\n")

    print(f"Written {output_file}")

    # Verify and report statistics
    z_values = []
    with open(output_file, 'r') as f:
        lines = f.readlines()
        # Skip first 3 header lines
        for line in lines[3:]:
            z_values.append(float(line.strip()))

    z_min = min(z_values)
    z_max = max(z_values)
    print(f"  Terrain stats: min={z_min:.3f}, max={z_max:.3f}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Gaussian terrain file in Askervein format",
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
