#!/usr/bin/env python3
"""MOST reference height on flat stretched meshes.

Runs MOST_Zref_Stretched.i (flat periodic column, dz0 = 10 m stretched by 1.1,
erf.most.zref unset) through the three MOST lookup paths a flat stretched mesh
can take, plus a uniform 10 m column with its bottom at z = 0 and at z = 300 m,
and checks the surface friction velocity.

  case                mesh / terrain                  MOST lookup           zref
  fitted_plane        StretchedDz, fitted (auto)      set_k_indices_T        5 m
  fitted_interp       VariableDz, StaticFittedMesh    set_norm_positions_T   5 m
  immersed_stretched  StretchedDz, ImmersedForcing    set_k_indices_N        5 m
  uniform_reference   ConstantDz, dz = 10 m           set_k_indices_N        5 m
  uniform_raised_base ConstantDz, prob_lo[2] = 300 m  set_k_indices_N        5 m
  raised_base_k_arr   the same, erf.most.k_arr_in = 0 set_k_indices_N        5 m
  raised_base_dz_odd  ConstantDz, dz = 1000/96 m, base 300 m  set_k_indices_N  5.21 m

The initial wind is sheared: sounding_most_zref ramps it linearly from 10 m/s at
the ground to 15 m/s at 50 m, constant above.  ERF interpolates the sounding
linearly to the cell centers, so every cell below 50 m, and a linear
interpolation between two of them, holds exactly the ramp value at its height.
At t = 0 the plane-averaged u* must therefore be kappa * U(zref) / ln(zref / z0)
with the wind taken AT the reference height: a lookup that reports one height
but samples the wind of another cell gives a different u*.  From u*(0) and the
reported zref the script also inverts the ramp for the height the wind actually
came from.  The stretched no-terrain column must follow the uniform 10 m column
over the run, since both use the 5 m first cell center, and the raised column must
follow it to the 6 digits of hist.dat: zref is a height above the ground.  Sounding
heights are absolute, like the mesh, so the raised columns are given the sounding
raised by 300 m (SOUNDING_SHIFT): the same wind at the same height above the ground.

Five more cases put a Witch of Agnesi ridge (prob.custom_terrain_type = WoA,
40 m high, half-width 100 m) under the same column, so that the first cell's
height above the ground differs from column to column, and run the default
through each terrain-fitted lookup (see HILL_CASES).  The reported range of
reference heights must be half of each column's first-cell thickness, from the
BTF mapping z = h + zeta (1 - h / z_top), and u*(0) the mean over the columns
of the log law with the wind of each column's first cell.  That wind is what
ERF puts there: the sounding is first resampled onto the flat level-0 cell
centres (ERF_InputSoundingData.H), each x face is then given the wind of that
table at its own physical height, and the first-cell wind of a column is the
mean of its two faces (erf_face_wind, column_winds).  On the ridge the crest
face lies between two resampled nodes that straddle the ramp's kink at 50 m, so
its wind is 1 % below the ramp value; the flat cases are exact with the ramp
only because their cell centres are the resampled nodes.  Before #4935 these
cases queried 10 m above the ground and took the wind of the second cell.

Before the fix the two fitted cases aborted at start-up (the default 10 m query
lies exactly on the top face of the first cell) and immersed_stretched reported
zref = 55.3 m, half of the uniform (prob_hi - prob_lo)/nz spacing, while taking
the wind of the first cell (5 m).  Before the reference height was made
wall-relative, the raised column reported zref = 305 m, the absolute height, and
gave u*(0) = 0.537 for 1.10.
Pass --old-exe to rerun a binary built before
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
    ("fitted_plane",       [],                                                      5.0, 15.5),
    ("fitted_interp",      ["erf.terrain_type=StaticFittedMesh"],                   5.0, 15.5),
    ("immersed_stretched", ["erf.terrain_type=ImmersedForcing",
                            "eb2.geom_type=all_regular"],                           5.0, 15.5),
    ("uniform_reference",  ["erf.grid_stretching_ratio=0",
                            "amr.n_cell=4", "4", "443",
                            "geometry.prob_extent=400", "400", "4430"],             5.0, 15.0),
    ("uniform_raised_base", ["erf.grid_stretching_ratio=0",
                             "amr.n_cell=4", "4", "443",
                             "geometry.prob_extent=400", "400", "4430",
                             "geometry.prob_lo=0", "0", "300"],                    5.0, 15.0),
    ("raised_base_k_arr",   ["erf.grid_stretching_ratio=0",
                             "amr.n_cell=4", "4", "443",
                             "geometry.prob_extent=400", "400", "4430",
                             "geometry.prob_lo=0", "0", "300",
                             "erf.most.k_arr_in=0"],                               5.0, 15.0),
    # a raised base where prob_lo + dz/2 - prob_lo rounds below dz/2 (dz = 1000/96 m)
    ("raised_base_dz_odd",  ["erf.grid_stretching_ratio=0",
                             "amr.n_cell=4", "4", "96",
                             "geometry.prob_extent=400", "400", "1000",
                             "geometry.prob_lo=0", "0", "300"],             500.0 / 96.0, 1500.0 / 96.0),
]

# A ridge under the column, so that the first cell's height above the ground varies by
# column: the default reference must follow it (erf-model/ERF#4935).  The ridge is kept
# below the top of the sounding's wind ramp (50 m), so that the sampled cell still matters.
HILL_HMAX = 40.0    # prob.hmax [m]
HILL_L = 100.0      # prob.L [m]
HILL_OPTS = ["erf.terrain_type=StaticFittedMesh", "prob.custom_terrain_type=WoA",
             "prob.hmax=%g" % HILL_HMAX, "prob.L=%g" % HILL_L]
# name, extra runtime options, MOST lookup exercised
HILL_CASES = [
    ("hill_default", [],                                             "set_norm_positions_T"),
    ("hill_plane",   ["erf.most.average_policy=0"],                  "set_k_indices_T"),
    ("hill_local",   ["erf.most.average_policy=1"],                  "set_k_indices_T"),
    ("hill_normal",  ["erf.most.average_policy=1",
                      "erf.most.use_normal_vector=true"],            "set_norm_indices_T"),
    ("hill_interp",  ["erf.most.average_policy=1",
                      "erf.most.use_interpolation=true"],            "set_z_positions_T"),
]
# the deck's grid, needed for the ridge's analytic first-cell heights
NX, DX = 4, 100.0           # amr.n_cell[0], prob_extent[0] / NX
NZ = 40                     # amr.n_cell[2]
DZ0 = 10.0                  # erf.initial_dz
RATIO = 1.1                 # erf.grid_stretching_ratio
ZTOP = 4425.925556817605    # prob_extent[2]: the top of the stretched column

# the raised columns against the uniform one: hist.dat holds 6 significant digits, so
# one unit in the last digit of u* ~ 1.1 is 9e-6
TOL_RAISED = 1.0e-5

# Sounding heights are absolute, like the mesh (the bottom of the domain is
# geometry.prob_lo[2], #4445), so a column whose base is raised gets the sounding raised
# with it: the same wind at the same height above its ground as the uniform column.
SOUNDING_SHIFT = {
    "uniform_raised_base": 300.0,
    "raised_base_k_arr":   300.0,
    "raised_base_dz_odd":  300.0,
}


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


def flat_levels():
    """Staggered heights of the flat level-0 column (erf.initial_dz stretched by the ratio)."""
    z, dz = [0.0], DZ0
    for _ in range(NZ):
        z.append(z[-1] + dz)
        dz *= RATIO
    return z


def erf_sounding_table(snd):
    """The wind table ERF interpolates from: the sounding resampled onto the flat cell
    centres, with the surface value first and the value at the top last."""
    zlev = flat_levels()
    z = [0.0] + [0.5 * (zlev[k] + zlev[k + 1]) for k in range(NZ)] + [zlev[-1]]
    return z, [wind_at(snd, zq) for zq in z]


def column_winds(snd):
    """First-cell wind of each ridge column as ERF initialises it: each x face takes the
    resampled table's wind at its own height, h_n + dz0/2 (1 - h_n / ztop), and the
    cell-centred wind is the mean of the column's two faces."""
    table = erf_sounding_table(snd)
    xc = 0.5 * NX * DX
    h = [HILL_HMAX / (1.0 + ((n * DX - xc) / HILL_L) ** 2) for n in range(NX + 1)]
    face = [wind_at(table, hn + 0.5 * DZ0 * (1.0 - hn / ZTOP)) for hn in h]
    return [0.5 * (face[i] + face[i + 1]) for i in range(NX)]


def hill_columns():
    """Per column of the ridge decks: (first-cell height above the ground, absolute height of
    the first cell centre, the same for the second cell).

    ERF puts the ridge on the nodes, h(x) = hmax / (1 + ((x - xc) / L)^2) with xc the domain
    centre, and maps the stretched levels with BTF, z(h, zeta) = h + zeta (1 - h / ztop).
    The MOST lookups average the four nodes of a column, and the mapping is linear in h, so
    a column's faces are those of the mean nodal height."""
    xc = 0.5 * NX * DX
    h = [HILL_HMAX / (1.0 + ((n * DX - xc) / HILL_L) ** 2) for n in range(NX + 1)]
    cols = []
    for i in range(NX):
        hbar = 0.5 * (h[i] + h[i + 1])
        zlev1 = hbar + DZ0 * (1.0 - hbar / ZTOP)
        zlev2 = hbar + (DZ0 + DZ0 * 1.1) * (1.0 - hbar / ZTOP)
        zref = 0.5 * (zlev1 - hbar)
        cols.append((zref, hbar + zref, 0.5 * (zlev1 + zlev2)))
    return cols


def write_sounding(src, dst, shift):
    """Copy an input sounding with every height raised by shift (the surface line is kept)."""
    with open(src) as fin, open(dst, "w") as fout:
        lines = fin.read().splitlines()
        fout.write(lines[0] + "\n")
        for line in lines[1:]:
            cols = line.split()
            if not cols:
                continue
            cols[0] = repr(float(cols[0]) + shift)
            fout.write(" ".join(cols) + "\n")


def run_case(exe, mpi_cmd, workdir, name, opts, src_dir):
    """Run one configuration; return (log text, u* series or None, ok)."""
    rundir = os.path.join(workdir, name)
    os.makedirs(rundir, exist_ok=True)
    with open(os.path.join(src_dir, DECK)) as fin, open(os.path.join(rundir, DECK), "w") as fout:
        fout.write(fin.read())
    write_sounding(os.path.join(src_dir, SOUNDING), os.path.join(rundir, SOUNDING),
                   SOUNDING_SHIFT.get(name, 0.0))
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


def reported_zref_range(log):
    """(smallest, largest) reference height from the log, or None."""
    m = re.search(r"MOST reference height at level 0: ([-+0-9.eE]+)(?: to ([-+0-9.eE]+))?", log)
    if m:
        lo = float(m.group(1))
        return (lo, float(m.group(2)) if m.group(2) else lo)
    return None


def reported_zref(log):
    """Reference height from the log, or None."""
    rng = reported_zref_range(log)
    if rng:
        return rng[0]
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

    # The ridge: the reference heights must span the columns' first-cell half-heights
    # and u*(0) must be the mean over the columns of the log law with each column's
    # first-cell wind (plane averaging takes the mean wind, local averaging each column's)
    cols = hill_columns()
    zref_lo, zref_hi = min(c[0] for c in cols), max(c[0] for c in cols)
    # the wind of the second cell must give a visibly different u* (ramp estimate)
    hill_right = sum(KAPPA * wind_at(snd, c[1]) / math.log(c[0] / Z0) for c in cols) / len(cols)
    hill_wrong = sum(KAPPA * wind_at(snd, c[2]) / math.log(c[0] / Z0) for c in cols) / len(cols)
    hill_sens = abs(hill_wrong - hill_right) / hill_right
    if hill_sens < MIN_SENSITIVITY:
        failures.append("ridge: sounding too uniform, the second cell's wind changes u* by only %.1e"
                        % hill_sens)
    # the first-cell wind of each column as ERF initialises it, and the expected u*(0):
    # plane averaging takes the mean wind, local averaging each column's own
    winds = column_winds(snd)
    umean = sum(winds) / len(winds)
    hill_expect = {}
    for name, _, _ in HILL_CASES:
        uw = [umean] * len(winds) if name == "hill_plane" else winds
        hill_expect[name] = sum(KAPPA * u / math.log(c[0] / Z0) for u, c in zip(uw, cols)) / len(cols)
    hill_new = {}
    for name, extra, _ in HILL_CASES:
        log, ustar, ok = run_case(args.exe, args.mpi_cmd, os.path.join(workdir, "new"),
                                  name, HILL_OPTS + extra, src_dir)
        rng = reported_zref_range(log)
        hill_new[name] = (ok, rng, ustar)
        if not ok or not ustar:
            failures.append("%s: run %s" % (name, abort_reason(log)))
            continue
        if (rng is None or abs(rng[0] - zref_lo) > TOL_ZREF * zref_lo
                or abs(rng[1] - zref_hi) > TOL_ZREF * zref_hi):
            failures.append("%s: zref %s, expected %.6g to %.6g" % (name, rng, zref_lo, zref_hi))
        expect = hill_expect[name]
        u0 = ustar[0][1]
        if abs(u0 - expect) > TOL_LOGLAW * expect:
            failures.append("%s: u*(0) %.6g, log law with each column's first-cell wind %.6g"
                            % (name, u0, expect))
        if any((not math.isfinite(u)) or u <= 0.0 for _, u in ustar):
            failures.append("%s: non-finite or non-positive u*" % name)

    old = {}
    hill_old = {}
    if args.old_exe:
        for name, opts, _, _ in CASES:
            log, ustar, ok = run_case(args.old_exe, args.mpi_cmd, os.path.join(workdir, "old"),
                                      name, opts, src_dir)
            old[name] = (ok, reported_zref(log), ustar, abort_reason(log))
        for name, extra, _ in HILL_CASES:
            log, ustar, ok = run_case(args.old_exe, args.mpi_cmd, os.path.join(workdir, "old"),
                                      name, HILL_OPTS + extra, src_dir)
            hill_old[name] = (ok, reported_zref_range(log), ustar, abort_reason(log))

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

    print()
    print("Witch of Agnesi ridge, %g m high, half-width %g m, under the same column:" % (HILL_HMAX, HILL_L))
    print("first-cell half heights %.5g to %.5g m; the ridge's second-cell wind changes u* by %.1e"
          % (zref_lo, zref_hi, hill_sens))
    header = "%-13s %-21s %18s %10s %10s %10s" % ("case", "lookup", "zref [m]", "u*(0)", "log law", "rel err")
    if hill_old:
        header += "   old run"
    print(header)
    for name, _, lookup in HILL_CASES:
        ok, rng, ustar = hill_new[name]
        expect = hill_expect[name]
        if ok and ustar and rng:
            u0 = ustar[0][1]
            row = "%-13s %-21s %8.5g to %7.5g %10.6g %10.6g %10.2e" % (
                name, lookup, rng[0], rng[1], u0, expect, abs(u0 - expect) / expect)
        else:
            row = "%-13s %-21s %18s %10s %10s %10s" % (name, lookup, "-", "failed", "-", "-")
        if hill_old:
            ook, orng, oustar, oreason = hill_old[name]
            if ook and oustar and orng:
                row += "   zref %.5g to %.5g, u*(0) %.6g" % (orng[0], orng[1], oustar[0][1])
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

    # The same column with its bottom at 300 m: the same u* at every step
    for raised in ("uniform_raised_base", "raised_base_k_arr"):
        ok_r, _, ur = new[raised]
        if not (ok_r and ok_u and ur and uu):
            continue
        worst = max(abs(a - b) / b for (_, a), (_, b) in zip(ur, uu))
        print()
        print("u* over the run, %s (bottom at 300 m) vs uniform_reference: largest rel diff %.2e (tol %.0e)"
              % (raised, worst, TOL_RAISED))
        if len(ur) != len(uu):
            failures.append("%s and uniform runs wrote %d and %d u* samples" % (raised, len(ur), len(uu)))
        if worst > TOL_RAISED:
            failures.append("%s u* departs from the uniform column by %.2e (tol %.0e)"
                            % (raised, worst, TOL_RAISED))

    print()
    if failures:
        for f in failures:
            print("FAIL: " + f)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
