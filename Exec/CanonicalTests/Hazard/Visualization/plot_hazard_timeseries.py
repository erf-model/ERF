#!/usr/bin/env python3
"""
plot_hazard_timeseries.py — ERF-Hazard diagnostic time series plots.

Reads plain-text CSV/DAT diagnostic files produced by ERF-Hazard
(dust_diag.dat, smoke_diag.dat, etc.) and produces publication-quality
line plots. No yt or AMReX dependency — only matplotlib and pandas.

Requirements:
    pip install matplotlib pandas numpy

Usage examples:

  # Terrain amplification: flat vs hill dust emission
  python plot_hazard_timeseries.py --mode terrain_amplification \
      --flat  path/to/HaboobFireFlat/dust_diag.dat \
      --hill  path/to/HaboobFireHill/dust_diag.dat \
      --pit   path/to/HaboobFirePit/dust_diag.dat

  # Single dust diagnostic file
  python plot_hazard_timeseries.py --mode dust_diag \
      --file path/to/dust_diag.dat --label "DustIntegration"

  # Smoke diagnostic
  python plot_hazard_timeseries.py --mode smoke_diag \
      --file path/to/smoke_diag.dat --label "HaboobFireHill"

  # Fire-dust coupling: all three interactions
  python plot_hazard_timeseries.py --mode coupling \
      --baseline path/to/FireDustBaseline/dust_diag.dat \
      --interaction1 path/to/FireDustInteraction1/dust_diag.dat \
      --interaction2 path/to/FireDustInteraction2/dust_diag.dat \
      --interaction3 path/to/FireDustInteraction3/dust_diag.dat

Output files are saved in the current directory.
"""

import argparse
import os
import sys

try:
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import numpy as np
except ImportError:
    sys.exit("ERROR: Install dependencies with: pip install matplotlib pandas numpy")

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "legend.fontsize": 11,
    "lines.linewidth": 2.0,
    "figure.dpi": 150,
})

COLORS = {
    "flat":         "#2196F3",   # blue
    "hill":         "#F44336",   # red
    "pit":          "#4CAF50",   # green
    "baseline":     "#9E9E9E",   # grey
    "interaction1": "#FF9800",   # orange — crust reduction
    "interaction2": "#9C27B0",   # purple — wind coupling
    "interaction3": "#E91E63",   # pink   — lofting
    "all":          "#212121",   # black  — all interactions
}

# ---------------------------------------------------------------------------
# File reader — handles both dust_diag.dat and smoke_diag.dat formats
# ---------------------------------------------------------------------------

def read_diag(path, label=""):
    """Read ERF diagnostic CSV/DAT file, skip comment lines starting with #."""
    if not os.path.isfile(path):
        sys.exit(f"ERROR: file not found: {path}")
    df = pd.read_csv(path, comment="#")
    df.columns = df.columns.str.strip()
    if label:
        df["_label"] = label
    print(f"  Loaded {path}: {len(df)} rows, columns: {list(df.columns)}")
    return df


# ---------------------------------------------------------------------------
# Plot modes
# ---------------------------------------------------------------------------

def plot_terrain_amplification(args):
    """
    Compare dust emission total over time for flat / hill / pit terrain.
    Shows how terrain geometry amplifies or suppresses emission.
    """
    datasets = {}
    if args.flat:
        datasets["Flat (reference)"] = (read_diag(args.flat, "flat"), COLORS["flat"])
    if args.hill:
        datasets["Hill +200 m"]      = (read_diag(args.hill, "hill"), COLORS["hill"])
    if args.pit:
        datasets["Pit -200 m"]       = (read_diag(args.pit,  "pit"),  COLORS["pit"])

    if not datasets:
        sys.exit("ERROR: provide at least one of --flat, --hill, --pit")

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    for label, (df, color) in datasets.items():
        t = df["time_s"]
        if "emission_total_kg_s" in df.columns:
            axes[0].plot(t, df["emission_total_kg_s"], color=color, label=label)
        if "ustar_max_m_s" in df.columns:
            axes[1].plot(t, df["ustar_max_m_s"], color=color, label=label)

    axes[0].set_ylabel("Dust emission total [kg/s]")
    axes[0].set_title("ERF-Hazard Phase 5 — Terrain Amplification of Dust Emission")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].set_ylabel("Max friction velocity u* [m/s]")
    axes[1].set_xlabel("Time [s]")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    outfile = "terrain_amplification.png"
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"Saved {outfile}")
    plt.close(fig)


def plot_dust_diag(args):
    """Plot all channels from a single dust_diag.dat file."""
    df = read_diag(args.file, args.label)
    t = df["time_s"]

    channels = [
        ("emission_total_kg_s",  "Emission total [kg/s]",    "#F44336"),
        ("deposition_total_kg_m2", "Deposition total [kg/m\u00b2]", "#2196F3"),
        ("ustar_max_m_s",        "Max u* [m/s]",             "#FF9800"),
        ("flux_max_kg_m2_s",     "Max emission flux [kg/m\u00b2/s]", "#9C27B0"),
        ("conc_sfc_max_kg_m3",   "Max sfc concentration [kg/m\u00b3]", "#4CAF50"),
    ]
    available = [(col, lbl, col_c) for col, lbl, col_c in channels if col in df.columns]

    if not available:
        sys.exit(f"No recognised columns found. Columns in file: {list(df.columns)}")

    n = len(available)
    fig, axes = plt.subplots(n, 1, figsize=(9, 2.5 * n), sharex=True)
    if n == 1:
        axes = [axes]

    for ax, (col, lbl, col_c) in zip(axes, available):
        ax.plot(t, df[col], color=col_c)
        ax.set_ylabel(lbl)
        ax.grid(True, alpha=0.3)

    axes[0].set_title(f"ERF-Hazard Dust Diagnostics — {args.label}")
    axes[-1].set_xlabel("Time [s]")
    fig.tight_layout()
    safe_label = args.label.replace(" ", "_")
    outfile = f"dust_diag_{safe_label}.png"
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"Saved {outfile}")
    plt.close(fig)


def plot_smoke_diag(args):
    """Plot smoke diagnostic channels."""
    df = read_diag(args.file, args.label)
    t = df["time_s"]

    channels = [
        ("smoke_src_max",       "Max smoke source [kg/m\u00b3/s]", "#FF5722"),
        ("smoke_conc_max_k0",   "Max sfc smoke [kg/m\u00b3]",      "#9C27B0"),
        ("smoke_total_mass",    "Total smoke mass [kg]",         "#3F51B5"),
    ]
    available = [(col, lbl, col_c) for col, lbl, col_c in channels if col in df.columns]

    if not available:
        # Fall back to whatever columns are present
        print(f"Standard smoke columns not found. Columns: {list(df.columns)}")
        available = [(col, col, "#333333") for col in df.columns if col != "time_s"]

    n = len(available)
    fig, axes = plt.subplots(n, 1, figsize=(9, 2.5 * n), sharex=True)
    if n == 1:
        axes = [axes]

    for ax, (col, lbl, col_c) in zip(axes, available):
        ax.plot(t, df[col], color=col_c)
        ax.set_ylabel(lbl)
        ax.grid(True, alpha=0.3)

    axes[0].set_title(f"ERF-Hazard Smoke Diagnostics — {args.label}")
    axes[-1].set_xlabel("Time [s]")
    fig.tight_layout()
    safe_label = args.label.replace(" ", "_")
    outfile = f"smoke_diag_{safe_label}.png"
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"Saved {outfile}")
    plt.close(fig)


def plot_coupling(args):
    """
    Compare dust emission across fire-dust coupling interactions.
    Shows incremental contribution of each interaction.
    """
    datasets = {}
    if args.baseline:
        datasets["Baseline (no coupling)"]   = (read_diag(args.baseline),     COLORS["baseline"])
    if args.interaction1:
        datasets["Interaction 1 (crust)"]     = (read_diag(args.interaction1), COLORS["interaction1"])
    if args.interaction2:
        datasets["Interaction 2 (wind u*)"]   = (read_diag(args.interaction2), COLORS["interaction2"])
    if args.interaction3:
        datasets["Interaction 3 (lofting)"]   = (read_diag(args.interaction3), COLORS["interaction3"])
    if args.all:
        datasets["All interactions"]           = (read_diag(args.all),          COLORS["all"])

    if not datasets:
        sys.exit("ERROR: provide at least one coupling dataset")

    fig, ax = plt.subplots(figsize=(9, 5))
    for label, (df, color) in datasets.items():
        if "emission_total_kg_s" in df.columns:
            ax.plot(df["time_s"], df["emission_total_kg_s"],
                    color=color, label=label)

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Dust emission total [kg/s]")
    ax.set_title("ERF-Hazard — Fire-Dust Coupling Interaction Contributions")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    outfile = "fire_dust_coupling.png"
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"Saved {outfile}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot ERF-Hazard diagnostic time series."
    )
    parser.add_argument(
        "--mode", required=True,
        choices=["terrain_amplification", "dust_diag", "smoke_diag", "coupling"],
        help="Plot mode"
    )

    # terrain_amplification
    parser.add_argument("--flat",  default=None, help="HaboobFireFlat dust_diag.dat")
    parser.add_argument("--hill",  default=None, help="HaboobFireHill dust_diag.dat")
    parser.add_argument("--pit",   default=None, help="HaboobFirePit  dust_diag.dat")

    # dust_diag / smoke_diag single file
    parser.add_argument("--file",  default=None, help="Single diagnostic file")
    parser.add_argument("--label", default="ERF-Hazard", help="Legend label")

    # coupling
    parser.add_argument("--baseline",     default=None)
    parser.add_argument("--interaction1", default=None)
    parser.add_argument("--interaction2", default=None)
    parser.add_argument("--interaction3", default=None)
    parser.add_argument("--all",          default=None,
                        help="All-interactions combined dust_diag.dat")

    args = parser.parse_args()

    if args.mode == "terrain_amplification":
        plot_terrain_amplification(args)
    elif args.mode == "dust_diag":
        if not args.file:
            sys.exit("ERROR: --file required for dust_diag mode")
        plot_dust_diag(args)
    elif args.mode == "smoke_diag":
        if not args.file:
            sys.exit("ERROR: --file required for smoke_diag mode")
        plot_smoke_diag(args)
    elif args.mode == "coupling":
        plot_coupling(args)

    print("Done.")


if __name__ == "__main__":
    main()
