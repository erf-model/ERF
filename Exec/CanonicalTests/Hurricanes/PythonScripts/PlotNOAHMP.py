import argparse
from pathlib import Path
import re
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate contour plots of TAU_EW and TAU_NS from NetCDF files."
    )
    parser.add_argument(
        "--lnd-dirs-path",
        type=str,
        default=".",
        help="Path to the directory containing lnd* folders (default: current directory)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Base directory containing the lnd* folders
    base_dir = Path(args.lnd_dirs_path).resolve()

    if not base_dir.exists():
        print(f"Error: Path '{base_dir}' does not exist.")
        return

    # Create target directories for output images inside the base directory
    out_tau_ew = base_dir / "TAU_EW"
    out_tau_ns = base_dir / "TAU_NS"

    out_tau_ew.mkdir(parents=True, exist_ok=True)
    out_tau_ns.mkdir(parents=True, exist_ok=True)

    # Find all matching directory paths (e.g., lnd10000, lnd11000)
    lnd_dirs = sorted(
        [d for d in base_dir.glob("lnd*") if d.is_dir()],
        key=lambda p: int(re.search(r"\d+", p.name).group())
        if re.search(r"\d+", p.name)
        else 0,
    )

    if not lnd_dirs:
        print(f"No directories matching 'lnd*' found in {base_dir}.")
        return

    print(f"Found {len(lnd_dirs)} directories to process in {base_dir}.")

    for lnd_dir in lnd_dirs:
        nc_path = lnd_dir / "Level_0.nc"

        if not nc_path.exists():
            print(f"Skipping {lnd_dir.name}: {nc_path.name} not found.")
            continue

        # Extract iteration number from directory name (e.g., 'lnd19000' -> '0019')
        iter_match = re.search(r"\d+", lnd_dir.name)
        iter_str = iter_match.group() if iter_match else lnd_dir.name
        iter_str = str(int(iter_str) // 1000).zfill(4)

        # Read NetCDF data
        with nc.Dataset(nc_path, mode="r") as ds:
            tau_ew = ds.variables["TAU_EW"][:]
            tau_ns = ds.variables["TAU_NS"][:]

        # Process and save TAU_EW
        plot_and_save(
            data=tau_ew,
            title=f"TAU_EW (iter={iter_str})",
            colorbar_label=r"$\tau_{EW}$",
            output_path=out_tau_ew / f"TAU_EW_iter_{iter_str}.png",
        )

        # Process and save TAU_NS
        plot_and_save(
            data=tau_ns,
            title=f"TAU_NS (iter={iter_str})",
            colorbar_label=r"$\tau_{NS}$",
            output_path=out_tau_ns / f"TAU_NS_iter_{iter_str}.png",
        )

        print(f"Processed: {lnd_dir.name} (iter={iter_str})")

    print("\nProcessing finished successfully!")


def plot_and_save(data, title, colorbar_label, output_path):
    fig, ax = plt.subplots(figsize=(10, 8), dpi=150)

    # Dynamic min/max per individual plot
    vmin, vmax = np.min(data), np.max(data)

    # Handle flat/constant arrays gracefully
    if vmin == vmax:
        vmin -= 1e-5
        vmax += 1e-5

    # 2D contour mesh plot
    cf = ax.contourf(
        data, levels=50, cmap="coolwarm", vmin=vmin, vmax=vmax, extend="both"
    )

    # Colorbar and Labels
    cbar = fig.colorbar(cf, ax=ax, orientation="vertical", pad=0.02)
    cbar.set_label(colorbar_label, fontsize=12)

    ax.set_title(f"iter={title.split('iter=')[1][:-1]}", fontsize=14, fontweight="bold")
    ax.set_xlabel("NX", fontsize=11)
    ax.set_ylabel("NY", fontsize=11)
    ax.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

