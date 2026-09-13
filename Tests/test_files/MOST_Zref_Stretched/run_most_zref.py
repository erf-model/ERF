#!/usr/bin/env python3
"""MOST reference height on flat stretched meshes.

Runs MOST_Zref_Stretched.i (flat periodic column, dz0 = 10 m stretched by 1.1,
erf.most.zref unset) through the three MOST lookup paths a flat stretched mesh
can take, plus a uniform 10 m column, and checks the surface friction velocity.

  case                mesh / terrain                  MOST lookup           zref
  fitted_plane        StretchedDz, fitted (auto)      set_k_indices_T       15.5 m
  fitted_interp       VariableDz, StaticFittedMesh    set_norm_positions_T  10 m
  immersed_stretched  StretchedDz, ImmersedForcing    set_k_indices_N        5 m
  uniform_reference   ConstantDz, dz = 10 m           set_k_indices_N        5 m

The initial wind is sheared: sounding_most_zref ramps it linearly from 10 m/s at
the ground to 15 m/s at 50 m, constant above.  ERF interpolates the sounding
linearly to the cell centers, so every cell below 50 m, and a linear
interpolation between two of them, holds exactly the ramp value at its height.
At t = 0 the plane-averaged u* must therefore be kappa * U(zref) / ln(zref / z0)
with the wind taken AT the reference height: a lookup that reports one height
but samples the wind of another cell gives a different u*.  From u*(0) and the
reported zref the script also inverts the ramp for the height the wind actually
came from.  The stretched no-terrain column must follow the uniform 10 m column
over the run, since both use the 5 m first cell center.

Before the fix the two fitted cases aborted at start-up (the default 10 m query
lies exactly on the top face of the first cell) and immersed_stretched reported
zref = 55.3 m, half of the uniform (prob_hi - prob_lo)/nz spacing, while taking
the wind of the first cell (5 m).  Pass --old-exe to rerun a binary built before
the fix and tabulate old against new.
"""

import argparse
import math
import os
import re
import shlex
import subprocess
import sys

KAPPA = 0.41      # erf.most von Karman constant (default)
Z0 = 0.1          # erf.most.z0 in the deck
DECK = "MOST_Zref_Stretched.i"
SOUNDING = "sounding_most_zref"

# relative tolerances: hist.dat holds 6 significant digits
TOL_ZREF = 1.0e-6
TOL_LOGLAW = 2.0e-5
TOL_UNIFORM = 1.0e-3
# the sampled wind must be told apart from the next cell's by far more than TOL_LOGLAW
MIN_SENSITIVITY = 50.0 * TOL_LOGLAW

# name, runtime options, expected zref [m], height of the neighbouring cell center [m]
CASES = [
    ("fitted_plane",       [],                                                     15.5, 27.05),
    ("fitted_interp",      ["erf.terrain_type=StaticFittedMesh"],                  10.0, 15.5),
    ("immersed_stretched", ["erf.terrain_type=ImmersedForcing",
                            "eb2.geom_type=all_regular"],                           5.0, 15.5),
    ("uniform_reference",  ["erf.grid_stretching_ratio=0",
                            "amr.n_cell=4", "4", "443",
                            "geometry.prob_extent=400", "400", "4430"],             5.0, 15.0),
]


def read_sounding(path):
    """Heights and wind speeds of an input sounding (first line is the surface)."""
    z, spd = [], []
    with open(path) as f:
        f.readline()
        for line in f:
            c = line.split()
            if len(c) >= 5:
                z.append(float(c[0]))
                spd.append(math.hypot(float(c[3]), float(c[4])))
    return z, spd


def wind_at(snd, zq):
    """Linear interpolation of the sounding wind speed, as ERF fills the cells."""
    z, spd = snd
    if zq <= z[0]:
        return spd[0]
    for k in range(len(z) - 1):
        if z[k] <= zq <= z[k + 1]:
            return spd[k] + (spd[k + 1] - spd[k]) * (zq - z[k]) / (z[k + 1] - z[k])
    return spd[-1]


def sampled_height(snd, u_star, zref):
    """Height whose sounding wind gives u_star with the reported zref (lowest match), or None."""
    target = u_star * math.log(zref / Z0) / KAPPA
    z, spd = snd
    for k in range(len(z) - 1):
        lo, hi = spd[k], spd[k + 1]
        if min(lo, hi) - 1e-9 <= target <= max(lo, hi) + 1e-9 and hi != lo:
            return z[k] + (target - lo) / (hi - lo) * (z[k + 1] - z[k])
    return None


def run_case(exe, mpi_cmd, workdir, name, opts, src_dir):
    """Run one configuration; return (log text, u* series or None, ok)."""
    rundir = os.path.join(workdir, name)
    os.makedirs(rundir, exist_ok=True)
    for f in (DECK, SOUNDING):
        with open(os.path.join(src_dir, f)) as fin, open(os.path.join(rundir, f), "w") as fout:
            fout.write(fin.read())
    hist = os.path.join(rundir, "hist.dat")
    if os.path.exists(hist):
        os.remove(hist)
    cmd = shlex.split(mpi_cmd) + [exe, DECK] + opts
    proc = subprocess.run(cmd, cwd=rundir, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True)
    with open(os.path.join(rundir, "run.log"), "w") as flog:
        flog.write(proc.stdout)
    ustar = None
    if os.path.exists(hist):
        ustar = []
        with open(hist) as fh:
            for line in fh:
                cols = line.split()
                if len(cols) >= 2 and not cols[0].isalpha():
                    try:
                        ustar.append((float(cols[0]), float(cols[1])))
                    except ValueError:
                        pass
    return proc.stdout, ustar, proc.returncode == 0


def reported_zref(log):
    """Reference height from the log, or None."""
    m = re.search(r"MOST reference height at level 0: ([-+0-9.eE]+)", log)
    if m:
        return float(m.group(1))
    m = re.search(r"Reference height for MOST set to ([-+0-9.eE]+)", log)
    if m:
        return float(m.group(1))
    return None


def abort_reason(log):
    m = re.search(r"Assertion `[^']*' failed, file \"[^\"]*/([^/\"]+)\", line (\d+), Msg: ([^\n]*)", log)
    if m:
        return "abort %s:%s %s" % (m.group(1), m.group(2), m.group(3).replace("!!!", "").strip()[:40])
    m = re.search(r"amrex::Abort::\d+::([^\n]*)", log)
    if m:
        return "abort %s" % m.group(1).replace("!!!", "").strip()[:60]
    return "failed"


def fmt_height(z):
    return "%.2f" % z if z is not None else "-"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--exe", required=True, help="erf_exec to test")
    parser.add_argument("--old-exe", default=None, help="erf_exec built before the fix (optional)")
    parser.add_argument("--mpi-cmd", default="", help="launcher prefix, e.g. 'mpiexec -n 1'")
    parser.add_argument("--workdir", default="most_zref_runs")
    args = parser.parse_args()

    src_dir = os.path.dirname(os.path.abspath(__file__))
    workdir = os.path.abspath(args.workdir)
    snd = read_sounding(os.path.join(src_dir, SOUNDING))
    loglaw = lambda zs, zr: KAPPA * wind_at(snd, zs) / math.log(zr / Z0)
    failures = []
    new = {}

    # The test only detects a wrong pairing if the neighbouring cell's wind differs
    for name, _, zexp, znext in CASES:
        sens = abs(loglaw(znext, zexp) - loglaw(zexp, zexp)) / loglaw(zexp, zexp)
        if sens < MIN_SENSITIVITY:
            failures.append("%s: sounding too uniform, wind %g m above zref changes u* by only %.1e"
                            % (name, znext, sens))

    for name, opts, zexp, _ in CASES:
        log, ustar, ok = run_case(args.exe, args.mpi_cmd, os.path.join(workdir, "new"),
                                  name, opts, src_dir)
        zref = reported_zref(log)
        new[name] = (ok, zref, ustar)
        if not ok or not ustar:
            failures.append("%s: run %s" % (name, abort_reason(log)))
            continue
        if zref is None or abs(zref - zexp) > TOL_ZREF * zexp:
            failures.append("%s: zref %s, expected %g" % (name, zref, zexp))
        u0 = ustar[0][1]
        expect = loglaw(zexp, zexp)
        if abs(u0 - expect) > TOL_LOGLAW * expect:
            zs = sampled_height(snd, u0, zref) if zref else None
            failures.append("%s: u*(0) %.6g, log law with the wind at %g m %.6g (wind taken at %s m)"
                            % (name, u0, zexp, expect, fmt_height(zs)))
        if any((not math.isfinite(u)) or u <= 0.0 for _, u in ustar):
            failures.append("%s: non-finite or non-positive u*" % name)

    old = {}
    if args.old_exe:
        for name, opts, _, _ in CASES:
            log, ustar, ok = run_case(args.old_exe, args.mpi_cmd, os.path.join(workdir, "old"),
                                      name, opts, src_dir)
            old[name] = (ok, reported_zref(log), ustar, abort_reason(log))

    print("MOST reference height on a flat stretched column (dz0 = 10 m, ratio 1.1)")
    print("initial wind %s; u* log law = %.2f * U(zref) / ln(zref / %.1f)"
          % (", ".join("%g m/s at %g m" % (s, z) for z, s in zip(*snd) if z <= 60.0), KAPPA, Z0))
    print("'wind at' is the height whose initial wind reproduces u*(0) with the reported zref")
    print()
    header = "%-19s %9s %10s %10s %10s %9s" % ("case", "zref [m]", "u*(0)", "log law", "rel err", "wind at")
    if old:
        header += "   old run"
    print(header)
    for name, _, zexp, _ in CASES:
        ok, zref, ustar = new[name]
        expect = loglaw(zexp, zexp)
        if ok and ustar:
            u0 = ustar[0][1]
            zs = sampled_height(snd, u0, zref) if zref else None
            row = "%-19s %9.4g %10.6g %10.6g %10.2e %9s" % (name, zref if zref is not None else float("nan"),
                                                          u0, expect, abs(u0 - expect) / expect, fmt_height(zs))
        else:
            row = "%-19s %9s %10s %10s %10s %9s" % (name, "-", "failed", "-", "-", "-")
        if old:
            ook, ozref, oustar, oreason = old[name]
            if ook and oustar:
                ozs = sampled_height(snd, oustar[0][1], ozref) if ozref else None
                row += "   zref %.4g, u*(0) %.6g, wind at %s m" % (ozref if ozref is not None else float("nan"),
                                                                  oustar[0][1], fmt_height(ozs))
            else:
                row += "   " + oreason
        print(row)

    # Stretched no-terrain column against the uniform 10 m column
    ok_s, _, us = new["immersed_stretched"]
    ok_u, _, uu = new["uniform_reference"]
    if ok_s and ok_u and us and uu:
        print()
        print("u* over the run: stretched (no terrain) vs uniform 10 m column")
        print("%8s %12s %12s %10s" % ("time", "stretched", "uniform", "rel diff"))
        worst = 0.0
        for (t, a), (_, b) in zip(us, uu):
            rel = abs(a - b) / b
            worst = max(worst, rel)
            print("%8.1f %12.6g %12.6g %10.2e" % (t, a, b, rel))
        if len(us) != len(uu):
            failures.append("stretched and uniform runs wrote %d and %d u* samples" % (len(us), len(uu)))
        if worst > TOL_UNIFORM:
            failures.append("stretched u* departs from the uniform column by %.2e (tol %.0e)"
                            % (worst, TOL_UNIFORM))

    print()
    if failures:
        for f in failures:
            print("FAIL: " + f)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
