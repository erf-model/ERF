#!/usr/bin/env python3
"""
Time the two-stream radiation solver against RRTMGP with everything else held
fixed, and plot the result.

Fairness: both solvers run from the same executable and share shared_settings,
so the only difference between the two runs at a given resolution is the
solver selected. Radiation is called every slow step in both, plotfiles and
per-solver logs are off, and the cost reported is the inclusive time AMReX's
TinyProfiler attributes to the solver's own region, not the wall time of the
run. Each configuration is repeated and the minimum is taken, since the
minimum is the least noisy estimator of a compute cost on a shared machine.

The cost is reported per step, not per call. The two-stream model is called
twice per step (pre_dycore runs the column sweep, post_dycore reuses it and
advances the surface balance), while RRTMGP is called once per step with
erf.rad_freq_in_steps = 1. A per-call figure would halve the two-stream cost.

Usage:
    ./run_timing_comparison.py --exe /path/to/erf_exec [--repeats 5]

The executable must be built with -DERF_ENABLE_TINY_PROFILE=ON, and with
-DERF_ENABLE_RRTMGP=ON for the RRTMGP half.
Without it the script still measures two-stream and says RRTMGP was skipped.
"""

import argparse
import csv
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

# Resolutions to sweep. Radiation cost scales with columns x levels, and the
# two solvers scale differently, so a single size would not be informative.
RESOLUTIONS = [(42, 4, 42), (84, 4, 42), (84, 4, 84), (168, 4, 84)]

# TinyProfiler regions at each solver's entry point. Both are read from the
# inclusive table, so time spent in nested profiled calls (for example AMReX
# communication) is charged to the solver that made them.
TWOSTREAM_REGION = "TwoStreamRadiation::advance()"
RRTMGP_REGION = "ERF::advance_radiation():RRTMGP"

# Two step counts; the difference isolates the marginal per-step cost.
STEPS_SHORT = 2
STEPS_LONG = 6


def parse_tiny_profiler(output, region):
    """Return (inclusive_seconds, ncalls) for a region, or None if absent.

    TinyProfiler prints an exclusive table and then an inclusive table; rows
    in both look like
      Name  NCalls  Min  Avg  Max  Max %
    The last matching row is the inclusive one. Avg is taken (equal to Min
    and Max on one rank).
    """
    found = None
    for line in output.splitlines():
        stripped = line.strip()
        # Per-region tables (BL_PROFILE_REGION) follow the main tables and
        # hold partial times; stop before them.
        if stripped.startswith("BEGIN REGION"):
            break
        # Timing rows end in a percentage; the memory table that follows
        # carries the same region names and is skipped. The name must be
        # followed by whitespace, so "<name>[sync]" rows do not match.
        if not stripped.startswith(region) or not stripped.endswith("%"):
            continue
        rest = stripped[len(region):]
        if not rest[:1].isspace():
            continue
        nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", rest)
        if len(nums) >= 4:
            found = (float(nums[2]), int(float(nums[0])))
    return found


def run_case(exe, case_dir, inputs, ncell, workdir, steps, extra):
    """Run one configuration and return the raw stdout."""
    for f in os.listdir(case_dir):
        src = os.path.join(case_dir, f)
        if os.path.isfile(src):
            shutil.copy(src, workdir)
    # print_threshold=0 keeps small regions out of the profiler's "Other" row.
    cmd = [exe, inputs, "amr.n_cell=%d %d %d" % ncell, "max_step=%d" % steps,
           "tiny_profiler.print_threshold=0"]
    cmd += extra
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-15:])
        raise RuntimeError("run failed (%s, %s):\n%s" % (inputs, ncell, tail))
    return proc.stdout


def measure_at(exe, case_dir, inputs, region, ncell, steps, repeats, extra):
    """Minimum inclusive time and call count over `repeats` runs of `steps` steps."""
    samples, ncalls = [], None
    for _ in range(repeats):
        with tempfile.TemporaryDirectory() as workdir:
            out = run_case(exe, case_dir, inputs, ncell, workdir, steps, extra)
        parsed = parse_tiny_profiler(out, region)
        if parsed is None:
            raise RuntimeError(
                "region %r not found in the profiler output for %s at %s; "
                "build ERF with -DERF_ENABLE_TINY_PROFILE=ON (the run sets "
                "tiny_profiler.print_threshold=0 so the region is not merged "
                "into the profiler's Other row)"
                % (region, inputs, ncell))
        samples.append(parsed[0])
        ncalls = parsed[1]
    return min(samples), ncalls


def measure(exe, case_dir, inputs, region, repeats, extra):
    """Return {ncell: (marginal_seconds_per_step, calls_long)} across the sweep.

    Each solver does one-time work on its first call: RRTMGP reads ~45 MB of
    k-distribution tables and both allocate scratch. Charging that to the
    per-step cost would flatter the solver with the larger setup. So each
    configuration is measured at two step counts and the reported figure is the
    marginal cost, (T_long - T_short) / (steps_long - steps_short), which
    cancels anything that happens once.
    """
    results = {}
    for ncell in RESOLUTIONS:
        t_short, c_short = measure_at(exe, case_dir, inputs, region, ncell,
                                      STEPS_SHORT, repeats, extra)
        t_long, c_long = measure_at(exe, case_dir, inputs, region, ncell,
                                    STEPS_LONG, repeats, extra)
        if c_long <= c_short:
            raise RuntimeError("call count did not increase with step count for %s" % inputs)
        marginal = (t_long - t_short) / (STEPS_LONG - STEPS_SHORT)
        results[ncell] = (marginal, c_long)
        print("    %-16s %-11s %8.3f ms/step (marginal, %d->%d calls)"
              % (inputs, "x".join(str(n) for n in ncell), marginal * 1e3,
                 c_short, c_long))
    return results


def write_csv(path, two, rr):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["nx", "ny", "nz", "columns", "cells",
                    "twostream_ms_per_step", "twostream_calls",
                    "rrtmgp_ms_per_step", "rrtmgp_calls", "cost_ratio"])
        for ncell in RESOLUTIONS:
            nx, ny, nz = ncell
            t, tc = two[ncell]
            if rr:
                r, rc = rr[ncell]
                speed = r / t if t > 0 else float("nan")
            else:
                r, rc, speed = "", "", ""
            w.writerow([nx, ny, nz, nx * ny, nx * ny * nz,
                        "%.4f" % (t * 1e3), tc,
                        ("%.4f" % (r * 1e3)) if r != "" else "", rc,
                        ("%.1f" % speed) if speed != "" else ""])


def make_plot(path, two, rr):
    """Plot what each solver costs and how that cost scales.

    The point of the figure is to let a user pick a solver for a given run, not
    to rank them: RRTMGP resolves hundreds of spectral g-points while the
    two-stream solver does one gray sweep, so the gap is the price of spectral
    detail rather than a defect. The right panel normalises by cell count to
    show that both scale linearly, which is the property that lets you
    extrapolate either one to your own problem size.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available; skipping the plot")
        return False

    cells = [nx * ny * nz for (nx, ny, nz) in RESOLUTIONS]
    tvals = [two[n][0] * 1e3 for n in RESOLUTIONS]          # ms per step
    tper = [two[n][0] / c * 1e6 for n, c in zip(RESOLUTIONS, cells)]  # us per cell

    ncols = 2
    fig, axes = plt.subplots(1, ncols, figsize=(11.5, 4.4))

    ax = axes[0]
    ax.plot(cells, tvals, "o-", color="#1f77b4", lw=2, ms=7, label="Two-stream (gray)")
    if rr:
        rvals = [rr[n][0] * 1e3 for n in RESOLUTIONS]
        ax.plot(cells, rvals, "s-", color="#d62728", lw=2, ms=7,
                label="RRTMGP (224 SW / 256 LW g-points)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("cells in the domain")
    ax.set_ylabel("radiation time per step [ms]")
    ax.set_title("Cost of one radiation update")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")

    ax2 = axes[1]
    ax2.plot(cells, tper, "o-", color="#1f77b4", lw=2, ms=7, label="Two-stream")
    if rr:
        rper = [rr[n][0] / c * 1e6 for n, c in zip(RESOLUTIONS, cells)]
        ax2.plot(cells, rper, "s-", color="#d62728", lw=2, ms=7, label="RRTMGP")
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("cells in the domain")
    ax2.set_ylabel("time per cell per step [us]")
    ax2.set_title("Cost per cell: both scale linearly")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)

    fig.suptitle("ERF radiation solvers: cost and scaling (single CPU core)",
                 fontsize=12, y=1.00)
    fig.text(0.5, -0.03,
             "Marginal cost per step, excluding one-time setup. Identical grid, "
             "timestep, sounding and moisture model; radiation called every step.\n"
             "The two solvers compute different things, so this is the cost of "
             "spectral detail, not a measure of quality. RRTMGP targets GPUs; "
             "this is CPU-only.",
             ha="center", fontsize=7.5, color="#444444")

    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print("  wrote %s" % path)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True, help="path to erf_exec")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--rrtmgp-data", default=None,
                    help="directory holding the four RRTMGP netCDF lookup tables")
    args = ap.parse_args()

    case_dir = os.path.dirname(os.path.abspath(__file__))
    exe = os.path.abspath(args.exe)
    if not os.path.isfile(exe):
        sys.exit("no executable at %s" % exe)

    print("Two-stream:")
    two = measure(exe, case_dir, "inputs_twostream", TWOSTREAM_REGION,
                  args.repeats, [])

    rr_extra = []
    if args.rrtmgp_data:
        rr_extra.append("erf.rrtmgp_file_path=%s" % os.path.abspath(args.rrtmgp_data))

    rr = None
    print("RRTMGP:")
    try:
        rr = measure(exe, case_dir, "inputs_rrtmgp", RRTMGP_REGION,
                     args.repeats, rr_extra)
    except RuntimeError as exc:
        print("  skipped: %s" % str(exc).splitlines()[0])
        print("  (an ERF built with -DERF_ENABLE_RRTMGP=ON is needed for this half)")

    csv_path = os.path.join(args.outdir, "radiation_timing_comparison.csv")
    write_csv(csv_path, two, rr)
    print("  wrote %s" % csv_path)
    make_plot(os.path.join(args.outdir, "radiation_timing_comparison.png"), two, rr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
