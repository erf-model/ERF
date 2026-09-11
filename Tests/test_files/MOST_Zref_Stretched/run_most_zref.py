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

The initial wind is a uniform 15 m/s, so at t = 0 the plane-averaged u* must be
the log-law value kappa * 15 / ln(zref / z0) at the reference height each path
reports.  The stretched no-terrain column must also follow the uniform 10 m
column over the run, since both use the 5 m first cell center.

Before the fix the two fitted cases aborted at start-up (the default 10 m query
lies exactly on the top face of the first cell) and immersed_stretched used
zref = 55.3 m, half of the uniform (prob_hi - prob_lo)/nz spacing, with the
first-cell wind (u* 0.974 instead of 1.572).  Pass --old-exe to rerun a binary
built before the fix and tabulate old against new.
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
U_INIT = 15.0     # uniform initial wind from the sounding [m/s]
DECK = "MOST_Zref_Stretched.i"
SOUNDING = "sounding_most_zref"

# relative tolerances: hist.dat holds 6 significant digits
TOL_ZREF = 1.0e-6
TOL_LOGLAW = 2.0e-5
TOL_UNIFORM = 1.0e-3

CASES = [
    ("fitted_plane",       [],                                                     15.5),
    ("fitted_interp",      ["erf.terrain_type=StaticFittedMesh"],                  10.0),
    ("immersed_stretched", ["erf.terrain_type=ImmersedForcing",
                            "eb2.geom_type=all_regular"],                           5.0),
    ("uniform_reference",  ["erf.grid_stretching_ratio=0",
                            "amr.n_cell=4", "4", "443",
                            "geometry.prob_extent=400", "400", "4430"],             5.0),
]


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
    return "failed"


def loglaw(zref):
    return KAPPA * U_INIT / math.log(zref / Z0)


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
    failures = []
    new = {}

    for name, opts, zexp in CASES:
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
        if abs(u0 - loglaw(zexp)) > TOL_LOGLAW * loglaw(zexp):
            failures.append("%s: u*(0) %.6g, log law at %g m %.6g" % (name, u0, zexp, loglaw(zexp)))
        if any((not math.isfinite(u)) or u <= 0.0 for _, u in ustar):
            failures.append("%s: non-finite or non-positive u*" % name)

    old = {}
    if args.old_exe:
        for name, opts, _ in CASES:
            log, ustar, ok = run_case(args.old_exe, args.mpi_cmd, os.path.join(workdir, "old"),
                                      name, opts, src_dir)
            old[name] = (ok, reported_zref(log), ustar, abort_reason(log))

    print("MOST reference height on a flat stretched column (dz0 = 10 m, ratio 1.1)")
    print("u* log law = %.2f * %.1f / ln(zref / %.1f)" % (KAPPA, U_INIT, Z0))
    print()
    header = "%-19s %9s %10s %10s %10s" % ("case", "zref [m]", "u*(0)", "log law", "rel err")
    if old:
        header += "   old run"
    print(header)
    for name, _, zexp in CASES:
        ok, zref, ustar = new[name]
        if ok and ustar:
            u0 = ustar[0][1]
            row = "%-19s %9.4g %10.6g %10.6g %10.2e" % (name, zref if zref is not None else float("nan"),
                                                      u0, loglaw(zexp), abs(u0 - loglaw(zexp)) / loglaw(zexp))
        else:
            row = "%-19s %9s %10s %10s %10s" % (name, "-", "failed", "-", "-")
        if old:
            ook, ozref, oustar, oreason = old[name]
            if ook and oustar:
                row += "   zref %.4g u*(0) %.6g" % (ozref if ozref is not None else float("nan"), oustar[0][1])
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
