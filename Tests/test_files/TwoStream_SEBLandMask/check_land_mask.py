#!/usr/bin/env python3
"""Check the land mask in the two-stream force-restore balance.

Two runs of TwoStream_SEBLandMask.i, one with erf.is_land = 1 and one with
erf.is_land = 0, each with a 2D plotfile per step. The deck is horizontally
uniform, so one slice in x carries the whole field. Asserts:

1. water run, every step: seb_t_sfc is --skin-temp and seb_q_sfc is
   --soil-water, exactly (a water column keeps its state);
2. land run: seb_t_sfc ends at least --min-change away from --skin-temp (the
   balance does evolve where it should);
3. step 1: seb_lw_net(land) - seb_lw_net(water) = sigma (T_sea^4 - T_skin^4)
   to --lw-rtol, where T_sea = t_surf (surf_pres / p0)^(R_d/c_p) from the water
   run's initial plotfile (the state the first sweep saw) and T_skin = --skin-temp. The deck's emissivity is one and
   both runs see the same sky at the first sweep, so nothing else differs;
4. the radiation diagnostics file has a finite T_s_mean on every post_dycore
   row of the land run and none on the water run (water columns are left out
   of the means);
5. with --nest-dir, a two-level run with water on level 0 and a level 1 of land
   over part of it: at the last step every level 0 column still holds
   --skin-temp and --soil-water (the average-down leaves a water column alone),
   while the level 1 skin has moved by at least --min-change.
"""

import argparse
import math
import os
import subprocess
import sys

P0 = 1.0e5                 # reference pressure [Pa]
RD = 287.0                 # gas constant of dry air [J/kg/K]
CP = 1004.5                # specific heat of dry air [J/kg/K]
SIGMA = 5.670374419e-8     # Stefan-Boltzmann constant [W/m^2/K^4]


def extract(fextract, plotfile, variable, out_path, fine_level=None):
    """Values of one variable along an x slice of a plotfile, from the finest level that
    covers each point, or from levels up to fine_level."""
    cmd = [fextract, '-d', '0', '-v', variable, '-s', out_path]
    if fine_level is not None:
        cmd += ['-f', str(fine_level)]
    cmd.append(plotfile)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"amrex_fextract failed on {plotfile} {variable} "
                           f"(exit {proc.returncode})\n{proc.stdout}{proc.stderr}")
    values = []
    with open(out_path) as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith('#'):
                values.append(float(line.split()[1]))
    if not values:
        raise RuntimeError(f"{plotfile}: no {variable} cells on the slice")
    return values


def skin_means(path):
    """T_s_mean of the post_dycore rows of a radiation diagnostics file."""
    with open(path) as handle:
        lines = [line.strip() for line in handle if line.strip()]
    header = [name.strip() for name in lines[0].split(',')]
    site, mean = header.index('call_site'), header.index('T_s_mean')
    out = []
    for line in lines[1:]:
        fields = [f.strip() for f in line.split(',')]
        if len(fields) == len(header) and fields[site] == 'post_dycore':
            out.append(float(fields[mean]))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fextract', required=True)
    parser.add_argument('--land-dir', required=True)
    parser.add_argument('--water-dir', required=True)
    parser.add_argument('--nest-dir', default=None)
    parser.add_argument('--steps', type=int, required=True)
    parser.add_argument('--prefix', default='plt2d')
    parser.add_argument('--digits', type=int, default=5)
    parser.add_argument('--skin-temp', type=float, default=299.0, help='erf.rad_t_sfc [K]')
    parser.add_argument('--soil-water', type=float, default=0.0095,
                        help='erf.radiation.seb_q_sfc_default')
    parser.add_argument('--min-change', type=float, default=1.0e-3, help='[K]')
    parser.add_argument('--min-lw-difference', type=float, default=5.0, help='[W/m^2]')
    parser.add_argument('--lw-rtol', type=float, default=1.0e-6)
    parser.add_argument('--sea-temp', type=float, default=297.1137,
                        help='301.5 (950/1000)^(R_d/c_p), from the first line of '
                             'input_sounding [K]')
    # 1e-4 is 0.03 K, or 33 Pa of surface pressure. The run is 5e-6 away: its surface
    # pressure is 1.6 Pa below 950 hPa (README). A cell-center pressure in place of the
    # surface pressure would be 8e-4 away.
    parser.add_argument('--sea-temp-rtol', type=float, default=1.0e-4)
    parser.add_argument('--start-rtol', type=float, default=1.0e-6,
                        help='plotfile 0 against the deck values')
    parser.add_argument('--diag-file', default='radiation_diag.dat')
    args = parser.parse_args()

    scratch = os.path.join(args.water_dir, 'land_mask_slices')
    os.makedirs(scratch, exist_ok=True)

    def field(directory, step, name):
        tag = os.path.basename(os.path.normpath(directory))
        plotfile = os.path.join(directory, f"{args.prefix}{step:0{args.digits}d}")
        return extract(args.fextract, plotfile, name, os.path.join(scratch, f"{tag}_{step}_{name}.dat"))

    failures = []

    # 1. The water columns keep their state: every plotfile holds the values of plotfile 0,
    # bit for bit, and plotfile 0 holds the deck's values (to rounding in the plotfile's
    # precision).
    for name, deck in (('seb_t_sfc', args.skin_temp), ('seb_q_sfc', args.soil_water)):
        start = field(args.water_dir, 0, name)
        if max(abs(v - deck) for v in start) > args.start_rtol * abs(deck):
            failures.append(f"water step 0: {name} is not the deck value {deck}")
        for step in range(1, args.steps + 1):
            now = field(args.water_dir, step, name)
            if len(now) != len(start) or any(a != b for a, b in zip(now, start)):
                failures.append(f"water step {step}: {name} is not what plotfile 0 holds")

    # 2. The land columns evolve.
    change = max(abs(v - args.skin_temp) for v in field(args.land_dir, args.steps, 'seb_t_sfc'))
    print(f"land: seb_t_sfc moved {change:.4e} K in {args.steps} steps")
    if change < args.min_change:
        failures.append(f"land: seb_t_sfc moved only {change:.3e} K; check 1 would be trivial")

    # 3. The longwave boundary over water is the sea surface temperature.
    # seb_lw_net of plotfile 1 comes from the sweep before step 1, which saw the initial
    # state: plotfile 0 holds that state's surface pressure.
    theta = field(args.water_dir, 0, 't_surf')
    pres = field(args.water_dir, 0, 'surf_pres')
    lw_land = field(args.land_dir, 1, 'seb_lw_net')
    lw_water = field(args.water_dir, 1, 'seb_lw_net')
    for i, (th, p, a, b) in enumerate(zip(theta, pres, lw_land, lw_water)):
        t_sea = th * (p / P0) ** (RD / CP)
        expected = SIGMA * (t_sea ** 4 - args.skin_temp ** 4)
        # A number from the sounding alone, not from the code's surface pressure.
        if abs(t_sea - args.sea_temp) > args.sea_temp_rtol * args.sea_temp:
            failures.append(f"column {i}: T_sea = {t_sea:.4f} K from the plotfile, "
                            f"{args.sea_temp} K from the sounding")
        if i == 0:
            print(f"step 1: T_sea = {t_sea:.4f} K, T_skin = {args.skin_temp} K; "
                  f"seb_lw_net land - water = {a - b:.6f} W/m^2, expected {expected:.6f}")
        if abs(expected) < args.min_lw_difference:
            failures.append(f"column {i}: the two boundaries differ by {expected:.3e} W/m^2, "
                            f"below {args.min_lw_difference}; the comparison would be trivial")
        if abs((a - b) - expected) > args.lw_rtol * abs(expected):
            failures.append(f"column {i}: seb_lw_net land - water = {a - b:.6f} W/m^2, "
                            f"expected sigma (T_sea^4 - T_skin^4) = {expected:.6f}")

    # 4. Water columns are left out of the means.
    land_means = skin_means(os.path.join(args.land_dir, args.diag_file))
    water_means = skin_means(os.path.join(args.water_dir, args.diag_file))
    if not land_means or not water_means:
        failures.append("no post_dycore rows in a radiation diagnostics file")
    if not all(math.isfinite(v) for v in land_means):
        failures.append("land: T_s_mean is not finite on some post_dycore row")
    if any(math.isfinite(v) for v in water_means):
        failures.append("water: T_s_mean is reported, though no column is in the balance")

    # 5. A coarse water column under fine land columns keeps its state.
    if args.nest_dir is not None:
        plotfile = os.path.join(args.nest_dir, f"{args.prefix}{args.steps:0{args.digits}d}")
        first = os.path.join(args.nest_dir, f"{args.prefix}{0:0{args.digits}d}")
        for name, deck in (('seb_t_sfc', args.skin_temp), ('seb_q_sfc', args.soil_water)):
            start = extract(args.fextract, first, name,
                            os.path.join(scratch, f"nest_coarse_start_{name}.slice"), fine_level=0)
            coarse = extract(args.fextract, plotfile, name,
                             os.path.join(scratch, f"nest_coarse_{name}.slice"), fine_level=0)
            if max(abs(v - deck) for v in start) > args.start_rtol * abs(deck):
                failures.append(f"nest step 0, level 0: {name} is not the deck value {deck}")
            if len(coarse) != len(start) or any(a != b for a, b in zip(coarse, start)):
                failures.append(f"nest, level 0 (water): {name} is not what plotfile 0 holds")
        both = extract(args.fextract, plotfile, 'seb_t_sfc',
                       os.path.join(scratch, 'nest_fine_seb_t_sfc.slice'))
        coarse = extract(args.fextract, plotfile, 'seb_t_sfc',
                         os.path.join(scratch, 'nest_coarse_seb_t_sfc_count.slice'), fine_level=0)
        if len(both) <= len(coarse):
            failures.append(f"nest: the slice holds {len(both)} points with level 1 and "
                            f"{len(coarse)} without; no level 1 column is on it")
        change = max(abs(v - args.skin_temp) for v in both)
        print(f"nest: level 1 seb_t_sfc moved {change:.4e} K in {args.steps} steps")
        if change < args.min_change:
            failures.append(f"nest: level 1 seb_t_sfc moved only {change:.3e} K; "
                            f"the average-down check would be trivial")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("PASS: water columns hold their state, take the sea surface temperature as the "
          "longwave boundary and are left out of the means")
    return 0


if __name__ == '__main__':
    sys.exit(main())
