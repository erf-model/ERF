#!/usr/bin/env python3
"""Check the coupling of the two-stream surface energy balance with the surface layer.

The prognostic surface energy balance advances the skin temperature with
C_s dT_s/dt = R_net - H - LE - G - (restore).  The surface layer puts its own
sensible (H) and latent (LE) heat fluxes into the air.  With
erf.radiation.seb_turbulent_flux_source = surface_layer the balance takes those
same fluxes, so the ground loses what the air gains; with = defaults it takes
the constants seb_hfx_default and seb_lh_default (0 here), the old behaviour.
With erf.radiation.seb_surface_layer_uses_skin = true the surface layer also
takes its surface temperature from the balance's skin, so the coupling runs
both ways.

Given the 2D plotfiles of one run in each mode, written every step, this
asserts:

1. surface_layer (one-way) and two-way runs, every step: the balance's H and
   LE (seb_hfx, seb_lh) equal the surface layer's sensible_heat_flux and
   latent_heat_flux to round-off, and the fluxes are not small (so the
   equality is not 0 == 0);
2. defaults run: seb_hfx and seb_lh are the constant defaults;
3. one-way and two-way runs: the skin ends lower than in the defaults run by
   the energy the fluxes removed, sum over steps of dt (H + LE) / C_s, to
   within --budget-rtol;
4. one-way run: the surface layer's t_surf stays at --most-surf-temp;
5. two-way run, steps 2..N: t_surf(n) = T_s(n-1) (p0 / p_sfc)^(R_d/c_p), the
   skin the balance ended step n-1 with as a potential temperature, where
   p_sfc is the surface pressure at the start of step n.  The 2D output
   surf_pres is the pressure at the centre of the lowest cell; p_sfc adds
   rho g dz/2 with rho from the ideal gas law at T_s, which leaves an error
   well under --theta-rtol.

With --multilevel the same checks run on every level of every plotfile, column
by column along a slice in x through the middle row of level 0 in y (on a finer
level, the first fine row inside that coarse row; a level that covers only part
of the slice reports its own columns, and one that covers none of it fails), and
check 3 is skipped:
under a finer level the coarse skin is the average of the fine one, so a coarse
column's skin no longer follows its own fluxes alone. Check 5 then also skips a
level at the step it is created, which has no skin of its own from the step before.
Without --multilevel the deck is horizontally uniform and one slice carries the
whole field.
"""

import argparse
import os
import subprocess
import sys

P0 = 1.0e5          # reference pressure [Pa]
RD = 287.0          # gas constant of dry air [J/kg/K]
CP = 1004.5         # specific heat of dry air [J/kg/K]
GRAV = 9.81         # [m/s^2]


class CheckError(Exception):
    """A condition that fails the check (reported as FAIL, not as a traceback)."""


def plotfile_geometry(plotfile):
    """(prob_lo_x, [dx of each level], [refinement ratio of each level]) from a native
    plotfile's Header."""
    with open(os.path.join(plotfile, 'Header')) as handle:
        lines = [line.strip() for line in handle]
    base = 2 + int(lines[1])          # version, ncomp, then the ncomp names
    finest = int(lines[base + 2])     # dim, time, finest level
    prob_lo_x = float(lines[base + 3].split()[0])
    ratios = [int(r) for r in lines[base + 5].split()] if finest > 0 else []
    dx = [float(lines[base + 8 + lev].split()[0]) for lev in range(finest + 1)]
    return prob_lo_x, dx, ratios


def extract(fextract, plotfile, variable, out_path, level=0):
    """[(x, value)] of the cells of one level along an x slice.

    amrex_fextract fixes the slice's transverse index on its coarse level and scales
    it by the refinement ratio only from there, so -c L -f L would read fine row
    j = jloc (the coarse middle-row index taken as a fine index), a different row on
    every level. Slice from level 0 instead, which reports the uncovered cells of the
    coarser levels as well, and keep the rows at this level's cell centres: with even
    refinement ratios no coarser cell centre is one.
    """
    prob_lo_x, dx, ratios = plotfile_geometry(plotfile)
    if any(r % 2 for r in ratios[:level]):
        raise CheckError(f"{plotfile}: refinement ratios {ratios[:level]} below level {level} "
                         f"include an odd one, so its cells cannot be told from coarser ones")
    cmd = [fextract, '-d', '0', '-v', variable, '-c', '0', '-f', str(level),
           '-e', '-s', out_path, plotfile]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"amrex_fextract failed on {plotfile} {variable} (exit {proc.returncode})\n"
            f"  command: {' '.join(cmd)}\n{proc.stdout}{proc.stderr}")
    values = []
    with open(out_path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{out_path}: expected 'x value', got {line!r}")
            values.append((float(fields[0]), float(fields[1])))
    if level > 0:
        def on_level(x):
            offset = (x - prob_lo_x) / dx[level] - 0.5
            return abs(offset - round(offset)) < 1.0e-6
        values = [(x, v) for x, v in values if on_level(x)]
    if not values:
        raise CheckError(f"{plotfile} level {level}: no {variable} cells on the slice")
    return values


def finest_level(plotfile):
    """The finest level recorded in a native plotfile's Header."""
    with open(os.path.join(plotfile, 'Header')) as handle:
        lines = [line.strip() for line in handle]
    ncomp = int(lines[1])
    return int(lines[2 + ncomp + 2])


def mean(values):
    return sum(values) / len(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fextract', required=True, help='amrex_fextract executable')
    parser.add_argument('--surface-layer-dir', required=True,
                        help='run directory with seb_turbulent_flux_source = surface_layer')
    parser.add_argument('--defaults-dir', required=True,
                        help='run directory with seb_turbulent_flux_source = defaults')
    parser.add_argument('--two-way-dir', default=None,
                        help='run directory with seb_surface_layer_uses_skin = true')
    parser.add_argument('--prefix', default='plt2d')
    parser.add_argument('--digits', type=int, default=5)
    parser.add_argument('--steps', type=int, required=True, help='last step; a plotfile per step')
    parser.add_argument('--dt', type=float, required=True, help='fixed time step [s]')
    parser.add_argument('--heat-capacity', type=float, required=True,
                        help='erf.radiation.seb_surface_heat_capacity [J/m^2/K]')
    parser.add_argument('--dz', type=float, default=50.0, help='height of the lowest cell [m]')
    parser.add_argument('--most-surf-temp', type=float, default=301.5,
                        help='erf.most.surf_temp of the deck [K]')
    parser.add_argument('--hfx-default', type=float, default=0.0)
    parser.add_argument('--lh-default', type=float, default=0.0)
    parser.add_argument('--min-flux', type=float, default=1.0,
                        help='smallest acceptable |H| and |LE| [W/m^2]')
    parser.add_argument('--rtol', type=float, default=1.0e-10,
                        help='relative tolerance of the flux equality')
    parser.add_argument('--budget-rtol', type=float, default=0.05)
    parser.add_argument('--theta-rtol', type=float, default=1.0e-4,
                        help='relative tolerance of the two-way surface temperature')
    parser.add_argument('--multilevel', action='store_true',
                        help='check every level column by column; skip the budget check')
    parser.add_argument('--min-spread', type=float, default=0.0,
                        help='with --multilevel, smallest acceptable max - min of H over '
                             'the coarse columns [W/m^2]: a uniform field would not tell '
                             'columns apart')
    args = parser.parse_args()

    scratch = os.path.join(args.surface_layer_dir, 'flux_source_slices')
    os.makedirs(scratch, exist_ok=True)

    def plotfile(directory, step):
        return os.path.join(directory, f"{args.prefix}{step:0{args.digits}d}")

    def levels(directory, step):
        return range(finest_level(plotfile(directory, step)) + 1) if args.multilevel else [0]

    def column_field(directory, step, name, level=0):
        """[(x, value)] of one variable on one level."""
        tag = os.path.basename(os.path.normpath(directory))
        return extract(args.fextract, plotfile(directory, step), name,
                       os.path.join(scratch, f"{tag}_{step}_{level}_{name}.dat"), level)

    def field(directory, step, name, level=0):
        return [v for _, v in column_field(directory, step, name, level)]

    failures = []

    def check_fluxes_and_budget(directory, label):
        """Assertions 1 and 3 for one leg; returns nothing, appends failures."""
        removed = 0.0  # sum of dt (H + LE) over the steps, per unit area [J/m^2]
        for step in range(1, args.steps + 1):
            for level in levels(directory, step):
                where = f"{label} step {step} level {level}"
                for seb_name, sl_name in (('seb_hfx', 'sensible_heat_flux'),
                                          ('seb_lh', 'latent_heat_flux')):
                    seb = field(directory, step, seb_name, level)
                    sl = field(directory, step, sl_name, level)
                    if len(seb) != len(sl):
                        failures.append(f"{where}: {seb_name} and {sl_name} have different lengths")
                        continue
                    scale = max(abs(v) for v in sl)
                    if scale < args.min_flux:
                        failures.append(f"{where}: |{sl_name}| peaks at {scale:.3e} W/m^2, below "
                                        f"{args.min_flux}; the comparison would be trivial")
                    err = max(abs(a - b) for a, b in zip(seb, sl))
                    if err > args.rtol * max(scale, 1.0):
                        failures.append(f"{where}: {seb_name} differs from {sl_name} by {err:.3e} "
                                        f"W/m^2 (tolerance {args.rtol * max(scale, 1.0):.3e})")
                if args.multilevel and level == 0 and args.min_spread > 0.0:
                    h = field(directory, step, 'sensible_heat_flux', 0)
                    if max(h) - min(h) < args.min_spread:
                        failures.append(f"{where}: H spans only {max(h) - min(h):.3e} W/m^2 over "
                                        f"the columns, below {args.min_spread}; the columns are "
                                        f"not told apart")
            removed += args.dt * (mean(field(directory, step, 'seb_hfx')) +
                                  mean(field(directory, step, 'seb_lh')))

        if args.multilevel:
            return

        t_leg = mean(field(directory, args.steps, 'seb_t_sfc'))
        t_def = mean(field(args.defaults_dir, args.steps, 'seb_t_sfc'))
        cooling = t_def - t_leg
        expected = removed / args.heat_capacity
        print(f"{label}: sum dt (H + LE) = {removed:.6e} J/m^2; expected skin cooling "
              f"{expected:.6e} K, measured {cooling:.6e} K (T_s {t_leg:.6f} K, "
              f"{t_def:.6f} K with defaults)")
        if not expected > 0.0:
            failures.append(f"{label}: the fluxes removed no energy ({removed:.3e} J/m^2)")
        elif abs(cooling - expected) > args.budget_rtol * expected:
            failures.append(f"{label}: skin cooling {cooling:.6e} K is not the {expected:.6e} K "
                            f"the fluxes removed (relative error "
                            f"{abs(cooling - expected) / expected:.3e}, tolerance "
                            f"{args.budget_rtol})")

    # 1 and 3, one-way
    check_fluxes_and_budget(args.surface_layer_dir, 'one-way')

    # 2: the defaults leg keeps the constants
    for step in range(1, args.steps + 1):
        for level in levels(args.defaults_dir, step):
            for name, default in (('seb_hfx', args.hfx_default), ('seb_lh', args.lh_default)):
                worst = max(abs(v - default) for v in field(args.defaults_dir, step, name, level))
                if worst > 1.0e-12 * max(abs(default), 1.0):
                    failures.append(f"defaults step {step} level {level}: {name} is not the "
                                    f"default {default} (off by {worst:.3e})")

    # 4: without the two-way option the surface layer keeps its own temperature
    for step in range(1, args.steps + 1):
        for level in levels(args.surface_layer_dir, step):
            worst = max(abs(v - args.most_surf_temp)
                        for v in field(args.surface_layer_dir, step, 't_surf', level))
            if worst > 1.0e-10 * args.most_surf_temp:
                failures.append(f"one-way step {step} level {level}: t_surf moved off "
                                f"erf.most.surf_temp = {args.most_surf_temp} (by {worst:.3e} K)")

    if args.two_way_dir:
        # 1 and 3, two-way
        check_fluxes_and_budget(args.two_way_dir, 'two-way')

        # 5: t_surf(n) is the skin of step n-1 as a potential temperature
        # Counted per level and step: a level whose columns all went unpaired must
        # fail on its own, not hide behind the comparisons of another level.
        kappa = RD / CP
        worst_rel = 0.0
        compared = {}
        for step in range(2, args.steps + 1):
            for level in levels(args.two_way_dir, step):
                if level not in levels(args.two_way_dir, step - 1):
                    continue  # created this step: no skin of its own from the step before
                t_surf = dict(column_field(args.two_way_dir, step, 't_surf', level))
                t_skin = dict(column_field(args.two_way_dir, step - 1, 'seb_t_sfc', level))
                p_cc = dict(column_field(args.two_way_dir, step - 1, 'surf_pres', level))
                paired = [x for x in sorted(t_surf) if x in t_skin and x in p_cc]
                if not paired:
                    failures.append(f"two-way step {step} level {level}: none of its "
                                    f"{len(t_surf)} columns has a skin from step {step - 1} "
                                    f"to compare with")
                compared[level] = compared.get(level, 0) + len(paired)
                for x in paired:
                    rho = p_cc[x] / (RD * t_skin[x])
                    p_sfc = p_cc[x] + rho * GRAV * 0.5 * args.dz
                    theta = t_skin[x] * (P0 / p_sfc) ** kappa
                    rel = abs(t_surf[x] - theta) / theta
                    worst_rel = max(worst_rel, rel)
                    if rel > args.theta_rtol:
                        failures.append(f"two-way step {step} level {level} x {x:g}: t_surf "
                                        f"{t_surf[x]:.6f} K is not the step {step - 1} skin "
                                        f"{t_skin[x]:.6f} K as potential temperature {theta:.6f} K "
                                        f"(relative error {rel:.3e}, tolerance {args.theta_rtol})")
        if not compared:
            failures.append("two-way: no level had a skin from the step before to compare with")
        per_level = ', '.join(f"level {lev}: {n}" for lev, n in sorted(compared.items()))
        print(f"two-way: t_surf follows the previous step's skin as potential temperature "
              f"in {sum(compared.values())} column comparisons ({per_level}), worst relative "
              f"error {worst_rel:.3e}")

    if failures:
        for message in failures:
            print(f"FAIL: {message}")
        return 1
    print(f"PASS: the balance removed the surface layer's fluxes at all {args.steps} steps"
          + (", and the surface layer took the balance's skin" if args.two_way_dir else ""))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except CheckError as error:
        print(f"FAIL: {error}")
        sys.exit(1)
