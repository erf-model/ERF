#!/usr/bin/env python3
"""
Analyze WDM6 output and compare with WSM6
"""

import numpy as np
import matplotlib.pyplot as plt
import yt

def load_plotfile(path):
    """Load ERF plotfile with yt"""
    return yt.load(path)

def analyze_number_concentrations(ds):
    """
    Analyze WDM6 number concentrations
    Returns statistics and checks bounds
    """
    # Get domain data
    ad = ds.all_data()

    # Extract number concentrations (if they exist)
    try:
        nc = ad['nc'].to('1/cm**3')  # Cloud droplet number
        nr = ad['nr'].to('1/L')      # Rain drop number
        nn = ad['nn'].to('1/cm**3')  # Aerosol number

        print("=== Number Concentration Statistics ===")
        print(f"Cloud droplet (nc):")
        print(f"  Min: {nc.min():.2e} cm^-3")
        print(f"  Max: {nc.max():.2e} cm^-3")
        print(f"  Mean (where nc>0): {nc[nc>0].mean():.2e} cm^-3")

        print(f"\nRain drop (nr):")
        print(f"  Min: {nr.min():.2e} L^-1")
        print(f"  Max: {nr.max():.2e} L^-1")
        print(f"  Mean (where nr>0): {nr[nr>0].mean() if np.any(nr>0) else 0:.2e} L^-1")

        print(f"\nAerosol (nn):")
        print(f"  Min: {nn.min():.2e} cm^-3")
        print(f"  Max: {nn.max():.2e} cm^-3")
        print(f"  Mean: {nn.mean():.2e} cm^-3")

        # Check bounds
        ncmin = 1e1  # m^-3 = 1e-5 cm^-3
        nrmin = 1e-2  # m^-3

        if np.any(nc[nc>0] < ncmin * 1e-6):
            print("\n⚠️  WARNING: nc below ncmin threshold!")
        else:
            print("\n✓ nc bounds OK")

        if np.any(nr[nr>0] < nrmin * 1e-6):
            print("⚠️  WARNING: nr below nrmin threshold!")
        else:
            print("✓ nr bounds OK")

        if np.any(nn < 0):
            print("⚠️  WARNING: Negative nn values!")
        else:
            print("✓ nn non-negative")

        return {'nc': nc, 'nr': nr, 'nn': nn}

    except KeyError as e:
        print(f"Number concentration field not found: {e}")
        return None

def analyze_droplet_size(ds):
    """
    Calculate effective radius from qc and nc
    """
    ad = ds.all_data()

    try:
        qc = ad['qc'].v  # kg/kg
        nc = ad['nc'].v  # #/m^3
        rho = ad['density'].v  # kg/m^3

        # Droplet effective radius assuming liquid water density
        # r_eff = (3 * qc * rho / (4 * pi * rho_water * nc))^(1/3)
        rho_water = 1000.0  # kg/m^3

        # Only where clouds exist
        cloud_mask = (qc > 1e-6) & (nc > 10)  # qc > 1 g/kg, nc > 10/m^3

        if np.any(cloud_mask):
            r_eff = ((3 * qc[cloud_mask] * rho[cloud_mask]) /
                     (4 * np.pi * rho_water * nc[cloud_mask]))**(1/3)

            # Convert to microns
            r_eff_um = r_eff * 1e6

            print("\n=== Cloud Droplet Effective Radius ===")
            print(f"Min: {r_eff_um.min():.2f} μm")
            print(f"Max: {r_eff_um.max():.2f} μm")
            print(f"Mean: {r_eff_um.mean():.2f} μm")

            if r_eff_um.min() < 5 or r_eff_um.max() > 50:
                print("⚠️  WARNING: Droplet sizes outside typical range (5-50 μm)")
            else:
                print("✓ Droplet sizes reasonable")

            return r_eff_um
        else:
            print("\n⚠️  No clouds found (qc and nc too small)")
            return None

    except KeyError as e:
        print(f"Field not found for droplet size: {e}")
        return None

def compare_wsm6_wdm6(wsm6_path, wdm6_path):
    """
    Compare WSM6 and WDM6 results
    """
    print("\n=== Comparing WSM6 vs WDM6 ===")

    ds_wsm6 = yt.load(wsm6_path)
    ds_wdm6 = yt.load(wdm6_path)

    ad_wsm6 = ds_wsm6.all_data()
    ad_wdm6 = ds_wdm6.all_data()

    # Compare moisture fields
    for var in ['qc', 'qr', 'qi', 'qs', 'qg']:
        try:
            wsm6_val = ad_wsm6[var].v
            wdm6_val = ad_wdm6[var].v

            diff_max = np.abs(wsm6_val - wdm6_val).max()
            rel_diff = diff_max / (wsm6_val.max() + 1e-20) * 100

            print(f"\n{var}:")
            print(f"  WSM6 max: {wsm6_val.max():.3e}")
            print(f"  WDM6 max: {wdm6_val.max():.3e}")
            print(f"  Max abs diff: {diff_max:.3e}")
            print(f"  Max rel diff: {rel_diff:.2f}%")

        except KeyError:
            print(f"  Field {var} not found")

    # Precipitation comparison
    try:
        # If rain accumulation fields exist
        rain_wsm6 = ad_wsm6['rain_accum'].v
        rain_wdm6 = ad_wdm6['rain_accum'].v

        print(f"\nTotal precipitation:")
        print(f"  WSM6: {rain_wsm6.sum():.3e}")
        print(f"  WDM6: {rain_wdm6.sum():.3e}")
        print(f"  Difference: {abs(rain_wsm6.sum() - rain_wdm6.sum())/rain_wsm6.sum()*100:.2f}%")

    except KeyError:
        print("  Precipitation accumulation fields not found")

def plot_vertical_profile(ds, output='wdm6_profile.png'):
    """
    Plot vertical profiles of WDM6 fields
    """
    ad = ds.all_data()

    # Get z-coordinate
    z = ad['z'].v

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    try:
        # Mixing ratios
        for i, var in enumerate(['qc', 'qr', 'qi']):
            ax = axes[0, i]
            val = ad[var].v
            ax.scatter(val, z, alpha=0.1, s=1)
            ax.set_xlabel(f'{var} (kg/kg)')
            ax.set_ylabel('Height (m)')
            ax.set_title(var.upper())
            ax.grid(True)

        # Number concentrations
        for i, (var, unit) in enumerate([('nc', 'cm$^{-3}$'),
                                          ('nr', 'L$^{-1}$'),
                                          ('nn', 'cm$^{-3}$')]):
            ax = axes[1, i]
            val = ad[var].v
            ax.scatter(val, z, alpha=0.1, s=1)
            ax.set_xlabel(f'{var} ({unit})')
            ax.set_ylabel('Height (m)')
            ax.set_title(var.upper())
            ax.set_xscale('log')
            ax.grid(True)

        plt.tight_layout()
        plt.savefig(output, dpi=150)
        print(f"\nPlot saved to {output}")

    except KeyError as e:
        print(f"Could not create profile plot: {e}")

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage:")
        print("  python analyze_wdm6.py <wdm6_plotfile>")
        print("  python analyze_wdm6.py <wdm6_plotfile> <wsm6_plotfile>")
        sys.exit(1)

    wdm6_path = sys.argv[1]

    print(f"Analyzing WDM6 output: {wdm6_path}")
    ds = load_plotfile(wdm6_path)

    # Analyze number concentrations
    nc_data = analyze_number_concentrations(ds)

    # Analyze droplet size
    r_eff = analyze_droplet_size(ds)

    # Create profile plot
    plot_vertical_profile(ds)

    # Compare with WSM6 if provided
    if len(sys.argv) >= 3:
        wsm6_path = sys.argv[2]
        compare_wsm6_wdm6(wsm6_path, wdm6_path)

    print("\n✓ Analysis complete")
