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

The deck is horizontally uniform, so one slice through the domain from
amrex_fextract carries the whole field.
"""

import argparse
import os
import subprocess
import sys

P0 = 1.0e5          # reference pressure [Pa]
RD = 287.0          # gas constant of dry air [J/kg/K]
CP = 1004.5         # specific heat of dry air [J/kg/K]
GRAV = 9.81         # [m/s^2]


def extract(fextract, plotfile, variable, out_path):
    """Run amrex_fextract along x on level 0 and return the values."""
    cmd = [fextract, '-d', '0', '-v', variable, '-c', '0', '-f', '0',
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
            values.append(float(fields[1]))
    if not values:
        raise ValueError(f"{out_path}: no data rows")
    return values


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
    args = parser.parse_args()

    scratch = os.path.join(args.surface_layer_dir, 'flux_source_slices')
    os.makedirs(scratch, exist_ok=True)

    def field(directory, step, name):
        plotfile = os.path.join(directory, f"{args.prefix}{step:0{args.digits}d}")
        tag = os.path.basename(os.path.normpath(directory))
        return extract(args.fextract, plotfile, name,
                       os.path.join(scratch, f"{tag}_{step}_{name}.dat"))

    failures = []

    def check_fluxes_and_budget(directory, label):
        """Assertions 1 and 3 for one leg; returns nothing, appends failures."""
        removed = 0.0  # sum of dt (H + LE) over the steps, per unit area [J/m^2]
        for step in range(1, args.steps + 1):
            for seb_name, sl_name in (('seb_hfx', 'sensible_heat_flux'),
                                      ('seb_lh', 'latent_heat_flux')):
                seb = field(directory, step, seb_name)
                sl = field(directory, step, sl_name)
                if len(seb) != len(sl):
                    failures.append(f"{label} step {step}: {seb_name} and {sl_name} "
                                    f"have different lengths")
                    continue
                scale = max(abs(v) for v in sl)
                if scale < args.min_flux:
                    failures.append(f"{label} step {step}: |{sl_name}| peaks at {scale:.3e} "
                                    f"W/m^2, below {args.min_flux}; the comparison would be "
                                    f"trivial")
                err = max(abs(a - b) for a, b in zip(seb, sl))
                if err > args.rtol * max(scale, 1.0):
                    failures.append(f"{label} step {step}: {seb_name} differs from {sl_name} "
                                    f"by {err:.3e} W/m^2 "
                                    f"(tolerance {args.rtol * max(scale, 1.0):.3e})")
            removed += args.dt * (mean(field(directory, step, 'seb_hfx')) +
                                  mean(field(directory, step, 'seb_lh')))

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
        for name, default in (('seb_hfx', args.hfx_default), ('seb_lh', args.lh_default)):
            worst = max(abs(v - default) for v in field(args.defaults_dir, step, name))
            if worst > 1.0e-12 * max(abs(default), 1.0):
                failures.append(f"defaults step {step}: {name} is not the default {default} "
                                f"(off by {worst:.3e})")

    # 4: without the two-way option the surface layer keeps its own temperature
    for step in range(1, args.steps + 1):
        worst = max(abs(v - args.most_surf_temp)
                    for v in field(args.surface_layer_dir, step, 't_surf'))
        if worst > 1.0e-10 * args.most_surf_temp:
            failures.append(f"one-way step {step}: t_surf moved off erf.most.surf_temp = "
                            f"{args.most_surf_temp} (by {worst:.3e} K)")

    if args.two_way_dir:
        # 1 and 3, two-way
        check_fluxes_and_budget(args.two_way_dir, 'two-way')

        # 5: t_surf(n) is the skin of step n-1 as a potential temperature
        kappa = RD / CP
        worst_rel = 0.0
        for step in range(2, args.steps + 1):
            t_surf = mean(field(args.two_way_dir, step, 't_surf'))
            t_skin = mean(field(args.two_way_dir, step - 1, 'seb_t_sfc'))
            p_cc = mean(field(args.two_way_dir, step - 1, 'surf_pres'))
            rho = p_cc / (RD * t_skin)
            p_sfc = p_cc + rho * GRAV * 0.5 * args.dz
            theta = t_skin * (P0 / p_sfc) ** kappa
            rel = abs(t_surf - theta) / theta
            worst_rel = max(worst_rel, rel)
            if rel > args.theta_rtol:
                failures.append(f"two-way step {step}: t_surf {t_surf:.6f} K is not the step "
                                f"{step - 1} skin {t_skin:.6f} K as potential temperature "
                                f"{theta:.6f} K (relative error {rel:.3e}, tolerance "
                                f"{args.theta_rtol})")
        print(f"two-way: t_surf follows the previous step's skin as potential temperature, "
              f"worst relative error {worst_rel:.3e}")

    if failures:
        for message in failures:
            print(f"FAIL: {message}")
        return 1
    print(f"PASS: the balance removed the surface layer's fluxes at all {args.steps} steps"
          + (", and the surface layer took the balance's skin" if args.two_way_dir else ""))
    return 0


if __name__ == '__main__':
    sys.exit(main())
