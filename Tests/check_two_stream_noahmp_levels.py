#!/usr/bin/env python3
"""Check two-stream radiation feeding Noah-MP on two levels (TwoStream_NoahMPLevels).

Level 0 is a grassland patch; level 1 covers its middle and, in the "own" and "nested"
legs, runs Noah-MP on a nested land file of its own that is barren, so its albedo and the
shortwave it absorbs differ from the grassland's. Both levels take one land step on the
first ERF step; plt2d00002 holds the result. Per leg:

own     level 1 spans the domain in z, sweeps its own columns, and runs Noah-MP itself.
        1. The run says level 1 runs the land model on its own setup file.
        2. On each level, the shortwave the sweep left absorbed at the ground (SW_surface,
           the step-1 row of the radiation CSV for that level) equals what that level's
           Noah-MP absorbed (sav + sag), within --sw-tol. On level 1 this ties the fine
           sweep, which reads the fine Noah-MP's albedo, to the fine Noah-MP, which
           integrates the fine sweep's forcing.
        3. Level 1's absorbed shortwave and albedo differ from level 0's, by more than
           --min-land-diff and 0.01: its land is its own, not level 0's interpolated.
nested  level 1 stops below the domain top, so it does not sweep.
        1. as above;
        2. the radiation CSV has no level-1 row (the level did not sweep);
        3. level 1's Noah-MP inputs (sw_flux_dn, lw_flux_dn, cos_zenith_angle) are those of
           level 0 in the parent cell, to round-off, and are not the missing-value fill;
        4. level 1's land differs from level 0's, as in "own".
interp  no level-1 land file.
        1. The run says level 1 takes its land state from level 0.
        2. Level 1's t_sfc, sav, sag and albedo equal level 0's in the parent cell, to
           round-off (the fields are uniform on level 0, so the interpolation is exact).

Fields are read along an x slice through the middle row of level 0 (on level 1, the first
fine row inside that coarse row) with amrex_fextract -c 0 -f L, keeping the rows at level
L's cell centres; with a refinement ratio of 2 no coarser centre is one of them.
"""

import argparse
import csv
import os
import subprocess
import sys


class CheckError(Exception):
    """A condition that fails the check (reported as FAIL, not as a traceback)."""


def plotfile_geometry(plotfile):
    """(prob_lo_x, [dx of each level], [refinement ratio of each level]) from the Header."""
    with open(os.path.join(plotfile, 'Header')) as handle:
        lines = [line.strip() for line in handle]
    base = 2 + int(lines[1])
    finest = int(lines[base + 2])
    prob_lo_x = float(lines[base + 3].split()[0])
    ratios = [int(r) for r in lines[base + 5].split()] if finest > 0 else []
    dx = [float(lines[base + 8 + lev].split()[0]) for lev in range(finest + 1)]
    return prob_lo_x, dx, ratios, finest


def extract(fextract, plotfile, variable, level, out_path):
    """{x: value} of the cells of one level on the slice."""
    prob_lo_x, dx, ratios, finest = plotfile_geometry(plotfile)
    if level > finest:
        raise CheckError(f"{plotfile} has no level {level} (finest level {finest})")
    if any(r % 2 for r in ratios[:level]):
        raise CheckError(f"{plotfile}: refinement ratios {ratios[:level]} include an odd one")
    cmd = [fextract, '-d', '0', '-v', variable, '-c', '0', '-f', str(level),
           '-e', '-s', out_path, plotfile]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise CheckError(f"amrex_fextract failed on {plotfile} {variable} (exit "
                         f"{proc.returncode}): {proc.stdout}{proc.stderr}")
    values = {}
    with open(out_path) as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].startswith('#'):
                continue
            x, v = float(fields[0]), float(fields[1])
            offset = (x - prob_lo_x) / dx[level] - 0.5
            if abs(offset - round(offset)) < 1.0e-6:
                values[x] = v
    if not values:
        raise CheckError(f"{plotfile} level {level}: no {variable} cells on the slice")
    return values


def parent_x(x, prob_lo_x, dx_coarse):
    """Cell centre of the level-0 cell that contains x."""
    i = int((x - prob_lo_x) // dx_coarse)
    return prob_lo_x + (i + 0.5) * dx_coarse


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--leg', required=True, choices=['own', 'nested', 'interp'])
    parser.add_argument('--run-dir', required=True)
    parser.add_argument('--fextract', required=True)
    parser.add_argument('--plotfile', default='plt2d00002')
    parser.add_argument('--csv', default='radiation_diag.csv')
    parser.add_argument('--sw-step', type=int, default=1,
                        help='radiation CSV step whose sweep the land step integrated on')
    parser.add_argument('--sw-tol', type=float, default=0.05,
                        help='absorbed shortwave agreement, sweep vs Noah-MP [W/m^2]')
    parser.add_argument('--min-land-diff', type=float, default=10.0,
                        help='smallest level-1 minus level-0 absorbed shortwave [W/m^2]')
    args = parser.parse_args()

    plotfile = os.path.join(args.run_dir, args.plotfile)
    scratch = os.path.join(args.run_dir, 'levels_slices')
    os.makedirs(scratch, exist_ok=True)
    prob_lo_x, dx, _, _ = plotfile_geometry(plotfile)
    failures = []

    def field(name, level):
        return extract(args.fextract, plotfile, name, level,
                       os.path.join(scratch, f"{name}_{level}.dat"))

    def absorbed(level):
        sav, sag = field('sav', level), field('sag', level)
        return {x: sav[x] + sag[x] for x in sav}

    def uniform_value(values, what):
        lo, hi = min(values.values()), max(values.values())
        if hi - lo > 1.0e-9 * max(abs(hi), 1.0):
            failures.append(f"{what} is not uniform on the slice ({lo} to {hi}); the "
                            f"comparison with the level's mean assumes it is")
        return sum(values.values()) / len(values)

    log_path = os.path.join(args.run_dir, 'simulation.log')
    if not os.path.exists(log_path):
        raise CheckError(f"{log_path} does not exist")
    with open(log_path) as handle:
        log = handle.read()
    own_line = "Noah-MP at level 1: runs the land model on its own setup file"
    interp_line = "Noah-MP at level 1: takes its land state from level 0"
    expected_line = interp_line if args.leg == 'interp' else own_line
    if expected_line not in log:
        failures.append(f"simulation.log does not say '{expected_line}'")

    rows = []
    csv_path = os.path.join(args.run_dir, args.csv)
    if os.path.exists(csv_path):
        with open(csv_path) as handle:
            rows = [r for r in csv.DictReader(handle)
                    if r['call_site'] == 'pre_dycore' and int(r['step']) == args.sw_step]
    sweep_sw = {int(r['level']): float(r['SW_surface']) for r in rows}

    if args.leg in ('own', 'nested'):
        land = {lev: uniform_value(absorbed(lev), f"level {lev} sav + sag") for lev in (0, 1)}
        albedo = {lev: uniform_value(field('albedo', lev), f"level {lev} albedo") for lev in (0, 1)}
        print(f"{args.leg}: Noah-MP absorbed shortwave level 0 {land[0]:.4f}, level 1 "
              f"{land[1]:.4f} W/m^2; albedo {albedo[0]:.4f}, {albedo[1]:.4f}")
        if abs(land[1] - land[0]) < args.min_land_diff or abs(albedo[1] - albedo[0]) < 0.01:
            failures.append(f"level 1's land is not its own: absorbed shortwave {land[1]:.4f} vs "
                            f"{land[0]:.4f} W/m^2, albedo {albedo[1]:.4f} vs {albedo[0]:.4f}")

    if args.leg == 'own':
        for lev in (0, 1):
            if lev not in sweep_sw:
                failures.append(f"no step-{args.sw_step} pre_dycore row for level {lev} in {args.csv}")
                continue
            diff = abs(sweep_sw[lev] - land[lev])
            print(f"own: level {lev} sweep SW_surface {sweep_sw[lev]:.4f} vs Noah-MP sav + sag "
                  f"{land[lev]:.4f} W/m^2 (difference {diff:.4f})")
            if diff > args.sw_tol:
                failures.append(f"level {lev}: the sweep left {sweep_sw[lev]:.4f} W/m^2 absorbed "
                                f"at the ground, Noah-MP absorbed {land[lev]:.4f} (difference "
                                f"{diff:.4f}, tolerance {args.sw_tol})")

    if args.leg == 'nested':
        if 1 in sweep_sw:
            failures.append(f"level 1 has a row in {args.csv}: it swept, so it is not a nested "
                            f"patch and this leg does not test the interpolated forcing")
        for name in ('sw_flux_dn', 'lw_flux_dn', 'cos_zenith_angle'):
            coarse, fine = field(name, 0), field(name, 1)
            worst = 0.0
            for x, v in fine.items():
                parent = coarse.get(parent_x(x, prob_lo_x, dx[0]))
                if parent is None:
                    failures.append(f"nested: no level-0 {name} in the parent cell of x {x:g}")
                    continue
                if not v > 0.0:
                    failures.append(f"nested: level-1 {name} at x {x:g} is {v}, not a forcing")
                worst = max(worst, abs(v - parent) / max(abs(parent), 1.0e-30))
            print(f"nested: level-1 {name} vs level 0 in the parent cell, worst relative "
                  f"difference {worst:.3e} over {len(fine)} cells")
            if worst > 1.0e-12:
                failures.append(f"nested: level-1 {name} is not level 0's in the parent cell "
                                f"(worst relative difference {worst:.3e})")

    if args.leg == 'interp':
        for name in ('t_sfc', 'sav', 'sag', 'albedo'):
            coarse, fine = field(name, 0), field(name, 1)
            worst = 0.0
            for x, v in fine.items():
                parent = coarse.get(parent_x(x, prob_lo_x, dx[0]))
                if parent is None:
                    failures.append(f"interp: no level-0 {name} in the parent cell of x {x:g}")
                    continue
                worst = max(worst, abs(v - parent))
            print(f"interp: level-1 {name} vs level 0 in the parent cell, worst difference "
                  f"{worst:.3e} over {len(fine)} cells")
            if worst > 1.0e-10 * max((abs(v) for v in coarse.values()), default=0.0):
                failures.append(f"interp: level-1 {name} is not level 0's (worst difference "
                                f"{worst:.3e})")

    if failures:
        for message in failures:
            print(f"FAIL: {message}")
        return 1
    print(f"PASS: {args.leg}")
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except CheckError as error:
        print(f"FAIL: {error}")
        sys.exit(1)
