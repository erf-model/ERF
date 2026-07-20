#!/usr/bin/env python3
"""
plot_hazard_fields.py — ERF-Hazard AMReX plotfile visualisation.

Reads ERF AMReX plotfile directories (plt_NNNNN) and produces PNG figures
for the Phase 5 canonical test cases: smoke plume plan view, vertical
cross-sections, wind speed maps, and dust emission maps.

Requirements:
    pip install yt matplotlib numpy

Usage:
    python plot_hazard_fields.py --plotdir path/to/plt_00020 --case HaboobFireHill
    python plot_hazard_fields.py --plotdir path/to/plt_00020 --case DustGaussianHill
    python plot_hazard_fields.py --plotdir path/to/plt_00020 --case HaboobFirePit
    python plot_hazard_fields.py --plotdir path/to/plt_00020 --case DustGaussianPit
    python plot_hazard_fields.py --plotdir path/to/plt_00020 --case HaboobFireFlat

Outputs (in current directory):
    <case>_smoke_plan.png       Smoke concentration plan view at k=0
    <case>_smoke_xz.png         Smoke vertical cross-section (x-z plane)
    <case>_wind_sfc.png         Horizontal wind speed at near-surface level
    <case>_dust_emission.png    Dust emission flux at surface (if available)
    <case>_wind_recirculation.png  Wind vectors in x-z (pit cases)

References:
    Koschmieder, H. (1924). Beitr. Phys. Atmos., 12, 33-55.
    yt project: https://yt-project.org/
"""

import argparse
import os
import sys

try:
    import yt
except ImportError:
    sys.exit("ERROR: yt not found. Install with: pip install yt")

try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import numpy as np
except ImportError:
    sys.exit("ERROR: matplotlib/numpy not found. Install with: pip install matplotlib numpy")


# ---------------------------------------------------------------------------
# Helper: safe slice plot — skips gracefully if field not in plotfile
# ---------------------------------------------------------------------------

def safe_slice(ds, axis, field, center, width, cmap, log_scale,
               cbar_label, title, outfile):
    """Create a yt SlicePlot and save to outfile. Returns True on success."""
    # Check field exists
    available = [f[1] for f in ds.field_list]
    if field not in available:
        print(f"  SKIP {outfile}: field '{field}' not in plotfile.")
        print(f"  Available fields: {available}")
        return False

    slc = yt.SlicePlot(ds, axis, ("boxlib", field),
                       center=center, width=width)
    slc.set_cmap(("boxlib", field), cmap)
    slc.set_log(("boxlib", field), log_scale)
    slc.set_colorbar_label(("boxlib", field), cbar_label)
    slc.set_title(("boxlib", field), title)
    slc.save(outfile)
    print(f"  Saved {outfile}")
    return True


# ---------------------------------------------------------------------------
# Per-field plot functions
# ---------------------------------------------------------------------------

def plot_smoke_plan(ds, case):
    """Plan view of smoke concentration at near-surface level."""
    dd = ds.domain_dimensions
    dom = ds.domain_right_edge - ds.domain_left_edge
    cx = float(ds.domain_left_edge[0]) + float(dom[0]) * 0.5
    cy = float(ds.domain_left_edge[1]) + float(dom[1]) * 0.5
    cz = float(ds.domain_left_edge[2]) + float(dom[2]) * 0.02  # near-surface

    safe_slice(
        ds, axis="z",
        field="smoke",
        center=[cx, cy, cz],
        width=[(float(dom[0]), "m"), (float(dom[1]), "m")],
        cmap="inferno",
        log_scale=False,
        cbar_label="Smoke density [kg/m\u00b3]",
        title=f"{case} — Smoke Plan View (near-surface)",
        outfile=f"{case}_smoke_plan.png"
    )


def plot_smoke_xz(ds, case):
    """Vertical cross-section of smoke in the x-z plane at domain y-centre."""
    dom = ds.domain_right_edge - ds.domain_left_edge
    cx = float(ds.domain_left_edge[0]) + float(dom[0]) * 0.5
    cy = float(ds.domain_left_edge[1]) + float(dom[1]) * 0.5
    cz = float(ds.domain_left_edge[2]) + float(dom[2]) * 0.5

    safe_slice(
        ds, axis="y",
        field="smoke",
        center=[cx, cy, cz],
        width=[(float(dom[0]), "m"), (float(dom[2]), "m")],
        cmap="inferno",
        log_scale=False,
        cbar_label="Smoke density [kg/m\u00b3]",
        title=f"{case} — Smoke x-z Cross-Section",
        outfile=f"{case}_smoke_xz.png"
    )


def plot_wind_sfc(ds, case):
    """Horizontal (x) wind speed at near-surface level."""
    dom = ds.domain_right_edge - ds.domain_left_edge
    cx = float(ds.domain_left_edge[0]) + float(dom[0]) * 0.5
    cy = float(ds.domain_left_edge[1]) + float(dom[1]) * 0.5
    cz = float(ds.domain_left_edge[2]) + float(dom[2]) * 0.02

    safe_slice(
        ds, axis="z",
        field="x_velocity",
        center=[cx, cy, cz],
        width=[(float(dom[0]), "m"), (float(dom[1]), "m")],
        cmap="RdBu_r",
        log_scale=False,
        cbar_label="U [m/s]",
        title=f"{case} — Surface Wind Speed",
        outfile=f"{case}_wind_sfc.png"
    )


def plot_wind_xz(ds, case):
    """Vertical cross-section of x-velocity — shows recirculation in pit cases."""
    dom = ds.domain_right_edge - ds.domain_left_edge
    cx = float(ds.domain_left_edge[0]) + float(dom[0]) * 0.5
    cy = float(ds.domain_left_edge[1]) + float(dom[1]) * 0.5
    cz = float(ds.domain_left_edge[2]) + float(dom[2]) * 0.5

    safe_slice(
        ds, axis="y",
        field="x_velocity",
        center=[cx, cy, cz],
        width=[(float(dom[0]), "m"), (float(dom[2]), "m")],
        cmap="RdBu_r",
        log_scale=False,
        cbar_label="U [m/s]",
        title=f"{case} — Wind x-z Cross-Section (recirculation)",
        outfile=f"{case}_wind_recirculation.png"
    )


def plot_dust_emission(ds, case):
    """Dust emission flux at surface."""
    dom = ds.domain_right_edge - ds.domain_left_edge
    cx = float(ds.domain_left_edge[0]) + float(dom[0]) * 0.5
    cy = float(ds.domain_left_edge[1]) + float(dom[1]) * 0.5
    cz = float(ds.domain_left_edge[2]) + float(dom[2]) * 0.01

    safe_slice(
        ds, axis="z",
        field="dust_emission_flux",
        center=[cx, cy, cz],
        width=[(float(dom[0]), "m"), (float(dom[1]), "m")],
        cmap="YlOrRd",
        log_scale=False,
        cbar_label="Dust emission flux [kg/m\u00b2/s]",
        title=f"{case} — Dust Emission Flux",
        outfile=f"{case}_dust_emission.png"
    )


def plot_dust_concentration(ds, case):
    """Total dust mass density in the atmosphere."""
    dom = ds.domain_right_edge - ds.domain_left_edge
    cx = float(ds.domain_left_edge[0]) + float(dom[0]) * 0.5
    cy = float(ds.domain_left_edge[1]) + float(dom[1]) * 0.5
    cz = float(ds.domain_left_edge[2]) + float(dom[2]) * 0.5

    safe_slice(
        ds, axis="y",
        field="rhoadv_dust",
        center=[cx, cy, cz],
        width=[(float(dom[0]), "m"), (float(dom[2]), "m")],
        cmap="YlOrRd",
        log_scale=False,
        cbar_label="Dust density [kg/m\u00b3]",
        title=f"{case} — Dust Concentration x-z",
        outfile=f"{case}_dust_xz.png"
    )


def plot_theta(ds, case):
    """Potential temperature vertical cross-section — shows thermal structure."""
    dom = ds.domain_right_edge - ds.domain_left_edge
    cx = float(ds.domain_left_edge[0]) + float(dom[0]) * 0.5
    cy = float(ds.domain_left_edge[1]) + float(dom[1]) * 0.5
    cz = float(ds.domain_left_edge[2]) + float(dom[2]) * 0.5

    safe_slice(
        ds, axis="y",
        field="theta",
        center=[cx, cy, cz],
        width=[(float(dom[0]), "m"), (float(dom[2]), "m")],
        cmap="coolwarm",
        log_scale=False,
        cbar_label="\u03b8 [K]",
        title=f"{case} — Potential Temperature x-z",
        outfile=f"{case}_theta_xz.png"
    )


# ---------------------------------------------------------------------------
# Case dispatch
# ---------------------------------------------------------------------------

CASES = [
    "HaboobFireHill",
    "HaboobFireFlat",
    "HaboobFirePit",
    "DustGaussianHill",
    "DustGaussianPit",
]


def make_plots(ds, case):
    print(f"\nGenerating plots for case: {case}")
    print(f"Fields in plotfile: {[f[1] for f in ds.field_list]}")

    if case in ("HaboobFireHill", "HaboobFireFlat"):
        plot_smoke_plan(ds, case)
        plot_smoke_xz(ds, case)
        plot_wind_sfc(ds, case)
        plot_dust_emission(ds, case)
        plot_dust_concentration(ds, case)
        plot_theta(ds, case)

    elif case == "HaboobFirePit":
        plot_wind_xz(ds, case)
        plot_smoke_xz(ds, case)
        plot_dust_concentration(ds, case)
        plot_theta(ds, case)

    elif case == "DustGaussianHill":
        plot_dust_emission(ds, case)
        plot_wind_sfc(ds, case)
        plot_dust_concentration(ds, case)

    elif case == "DustGaussianPit":
        plot_wind_xz(ds, case)
        plot_dust_emission(ds, case)
        plot_dust_concentration(ds, case)

    else:
        print(f"Unknown case '{case}' — plotting all available fields.")
        plot_smoke_plan(ds, case)
        plot_smoke_xz(ds, case)
        plot_wind_sfc(ds, case)
        plot_dust_emission(ds, case)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Visualise ERF-Hazard AMReX plotfile output."
    )
    parser.add_argument(
        "--plotdir", required=True,
        help="Path to AMReX plotfile directory (e.g. plt_00020)"
    )
    parser.add_argument(
        "--case", default="HaboobFireHill",
        choices=CASES,
        help="Test case name — controls which figures are generated"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.plotdir):
        sys.exit(f"ERROR: plotfile directory not found: {args.plotdir}")

    print(f"Loading plotfile: {args.plotdir}")
    ds = yt.load(args.plotdir)
    make_plots(ds, args.case)
    print("\nDone.")


if __name__ == "__main__":
    main()
