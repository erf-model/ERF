#!/usr/bin/env python3
"""Check the moisture coupling of the two-stream surface energy balance with the surface layer.

With erf.radiation.seb_surface_layer_uses_moisture the surface layer takes its land surface
mixing ratio as beta * q_sat(T_s, p_s) + (1 - beta) * q_air, with
beta = clamp((q_s - wilt) / (fc - wilt), 0, 1) from the balance's soil water q_s, and the
balance drains q_s by the latent heat flux it removes. With a soil type, beta comes instead
from Noah-MP's bare-soil resistance (and a canopy's with a vegetation type) in series with the
aerodynamic one, the bare soil evaporating from its pore air (Noah-MP's relative humidity
times q_sat). Given the 2D plotfiles (one per step) of eight runs this asserts the
following. The runs are: three on the linear factor alone whose soil starts and restores at
the wilting point (dry), at field capacity (wet) and in between (mid); mid's soil as a Noah-MP
soil type (bare); bare with a vegetation type at fraction 0 (bare_f0); bare at the wilting
point (bare_dry); bare under vegetation (veg); and a one-step veg without erf.most.z0
(tables).

1. dry: beta = 0, so the latent heat flux is zero at every step, while
2. wet: the latent heat flux is not small, so 1 is not 0 == 0;
3. step 1, when the three runs differ only in beta (same air, same skin):
   q_surf(wet) - q_surf(dry) is not small (the saturated surface is moister than the air),
   q_surf(mid) - q_surf(dry) = beta(mid) * (q_surf(wet) - q_surf(dry)) to round-off -- the
   linear blend with the dry run's q_surf as the air's -- and LE(mid) lies between LE(dry)
   and LE(wet). (LE itself is not exactly beta times the wet value: the moisture flux enters
   the surface layer's stability through the virtual heat flux.)
4. mid, every step n >= 2: the balance's water content follows its budget,
   q_s(n) = q_s(n-1) - dt * LE(n) / (L_v rho_w d_s) - dt * (q_s(n-1) - q_deep) / tau_q,
   with LE(n) the latent heat flux it removed (seb_lh): the water the air gains leaves the
   soil; and the soil has dried from step 1 to the last;
5. veg (bare's soil under Noah-MP's grassland, seb_vegetation_type): the canopy and soil
   resistances lower LE below mid's, which has the soil-water factor alone, but not to
   zero, at every step; and the water budget and drying of 4 hold.
6. bare: the soil resistance lowers LE below mid's linear factor at the same water content,
   but not to zero, at every step; and the water budget and drying of 4 hold;
7. bare_f0 = bare: a vegetation type at fraction 0 gives bare soil's q_surf and LE at every
   step, to round-off (the model is continuous as the vegetated fraction goes to 0);
8. job_info records the roughness each run used: the tables' value in the tables run
   (Noah-MP's Z0MVT and Z0SOIL by the vegetated fraction), the deck's erf.most.z0 in veg,
   each as the only erf.most.z0 entry;
9. bare_dry: at the wilting point the pore air is nearly dry, so the bare soil does not
   evaporate: LE <= 0 at every step (Noah-MP gives about -4 W/m^2 there).
"""

import argparse
import os
import subprocess
import sys

L_V = 2.5e6       # latent heat of vaporization [J/kg] (ERF_Constants.H)
RHO_W = 1000.0    # density of water [kg/m^3] (rhor, ERF_MicrophysicsConstants.H)


class CheckError(Exception):
    """A condition that fails the check (reported as FAIL, not as a traceback)."""


def mean_value(fextract, plotfile, variable, out_path):
    cmd = [fextract, '-d', '0', '-v', variable, '-e', '-s', out_path, plotfile]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise CheckError(f"amrex_fextract failed on {plotfile} {variable} "
                         f"(exit {proc.returncode}): {proc.stdout}{proc.stderr}")
    values = []
    with open(out_path) as handle:
        for line in handle:
            fields = line.split()
            if fields and not fields[0].startswith('#'):
                values.append(float(fields[1]))
    if not values:
        raise CheckError(f"{plotfile}: no {variable} cells on the slice")
    return sum(values) / len(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fextract', required=True)
    parser.add_argument('--dry-dir', required=True)
    parser.add_argument('--wet-dir', required=True)
    parser.add_argument('--mid-dir', required=True)
    parser.add_argument('--veg-dir', required=True)
    parser.add_argument('--bare-dir', required=True, help='soil type, no vegetation')
    parser.add_argument('--bare-f0-dir', required=True,
                        help='soil type and a vegetation type at fraction 0')
    parser.add_argument('--bare-dry-dir', required=True,
                        help='soil type at the wilting point, no vegetation')
    parser.add_argument('--continuity-rtol', type=float, default=1.0e-12,
                        help='relative tolerance of bare_f0 = bare')
    parser.add_argument('--tables-dir', required=True,
                        help='run without erf.most.z0 (roughness from the tables)')
    parser.add_argument('--tables-z0', type=float, required=True,
                        help="the tables' land roughness for the tables run [m]")
    parser.add_argument('--given-z0', type=float, required=True,
                        help="the deck's erf.most.z0 [m], which the veg run keeps")
    parser.add_argument('--mid-q', type=float, required=True, help='q_s and q_deep of the mid run')
    parser.add_argument('--wilt', type=float, required=True)
    parser.add_argument('--fc', type=float, required=True)
    parser.add_argument('--steps', type=int, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--depth', type=float, default=0.1, help='seb_moisture_layer_depth_m')
    parser.add_argument('--tau-q', type=float, default=86400.0,
                        help='seb_moisture_restore_timescale_s')
    parser.add_argument('--min-le', type=float, default=1.0,
                        help='smallest acceptable LE of the wet run [W/m^2]')
    parser.add_argument('--mix-rtol', type=float, default=1.0e-10,
                        help='relative tolerance of the surface mixing ratio relation')
    parser.add_argument('--budget-atol', type=float, default=1.0e-12,
                        help='water-content budget tolerance per step [m^3/m^3]')
    args = parser.parse_args()

    failures = []

    def series(directory, variable):
        scratch = os.path.join(directory, 'moisture_slices')
        os.makedirs(scratch, exist_ok=True)
        return [mean_value(args.fextract, os.path.join(directory, f"plt2d{step:05d}"), variable,
                           os.path.join(scratch, f"{variable}_{step}.dat"))
                for step in range(0, args.steps + 1)]

    # 1 and 2
    le_dry = series(args.dry_dir, 'latent_heat_flux')
    le_wet = series(args.wet_dir, 'latent_heat_flux')
    worst_dry = max(abs(v) for v in le_dry[1:])
    print(f"dry: largest |LE| over the steps {worst_dry:.3e} W/m^2")
    if worst_dry > 1.0e-6:
        failures.append(f"dry: LE reaches {worst_dry:.3e} W/m^2; at the wilting point the "
                        f"surface mixing ratio is the air's and LE must be zero")
    print(f"wet: LE at step 1 {le_wet[1]:.4f} W/m^2")
    if not le_wet[1] > args.min_le:
        failures.append(f"wet: LE at step 1 is {le_wet[1]:.3e} W/m^2, not above {args.min_le}; "
                        f"the dry check would be trivial")

    # 3
    beta_mid = min(1.0, max(0.0, (args.mid_q - args.wilt) / (args.fc - args.wilt)))
    q_dry = series(args.dry_dir, 'q_surf')[1]
    q_wet = series(args.wet_dir, 'q_surf')[1]
    q_mid = series(args.mid_dir, 'q_surf')[1]
    print(f"step 1 q_surf: dry {q_dry:.8f} (the air's), mid {q_mid:.8f}, wet {q_wet:.8f} kg/kg; "
          f"beta(mid) = {beta_mid:.4f}")
    if not q_wet - q_dry > 1.0e-4:
        failures.append(f"q_surf(wet) - q_surf(dry) is {q_wet - q_dry:.3e} kg/kg; the "
                        f"saturated surface must be moister than the air for check 3 to mean "
                        f"anything")
    else:
        expected = beta_mid * (q_wet - q_dry)
        error = abs((q_mid - q_dry) - expected)
        print(f"mid: q_surf(mid) - q_surf(dry) = {q_mid - q_dry:.8e}, beta * (wet - dry) = "
              f"{expected:.8e}")
        if error > args.mix_rtol * abs(expected):
            failures.append(f"mid: q_surf(mid) - q_surf(dry) = {q_mid - q_dry:.8e} is not "
                            f"beta * (q_surf(wet) - q_surf(dry)) = {expected:.8e}")
    le_mid = series(args.mid_dir, 'latent_heat_flux')
    print(f"step 1 LE: dry {le_dry[1]:.4f}, mid {le_mid[1]:.4f}, wet {le_wet[1]:.4f} W/m^2")
    if not le_dry[1] < le_mid[1] < le_wet[1]:
        failures.append("step 1: LE(mid) does not lie between LE(dry) and LE(wet)")

    # 4 (and the budget of 5)
    def check_budget(directory, label, q_deep, expect_drying=True):
        q_s = series(directory, 'seb_q_sfc')
        seb_lh = series(directory, 'seb_lh')
        worst = 0.0
        for n in range(2, args.steps + 1):
            expected = (q_s[n - 1] - args.dt * seb_lh[n] / (L_V * RHO_W * args.depth)
                        - args.dt * (q_s[n - 1] - q_deep) / args.tau_q)
            worst = max(worst, abs(q_s[n] - expected))
        print(f"{label}: soil water {q_s[1]:.8f} -> {q_s[args.steps]:.8f} m^3/m^3 over steps "
              f"1-{args.steps}; worst budget error {worst:.3e}")
        if worst > args.budget_atol:
            failures.append(f"{label}: the soil water does not follow its budget (worst error "
                            f"{worst:.3e}, tolerance {args.budget_atol})")
        if expect_drying and not q_s[1] - q_s[args.steps] > 0.0:
            failures.append(f"{label}: the soil did not dry ({q_s[1]} -> {q_s[args.steps]})")

    check_budget(args.mid_dir, 'mid', args.mid_q)

    # 5
    le_veg = series(args.veg_dir, 'latent_heat_flux')
    print(f"veg: LE at step 1 {le_veg[1]:.4f} W/m^2 (mid {le_mid[1]:.4f}); ratio over the steps "
          f"{min(v / m for v, m in zip(le_veg[1:], le_mid[1:])):.4f}-"
          f"{max(v / m for v, m in zip(le_veg[1:], le_mid[1:])):.4f}")
    for n in range(1, args.steps + 1):
        if not 0.0 < le_veg[n] < le_mid[n]:
            failures.append(f"veg step {n}: LE {le_veg[n]:.4f} W/m^2 is not between 0 and mid's "
                            f"{le_mid[n]:.4f}")
            break
    check_budget(args.veg_dir, 'veg', args.mid_q)

    # 6
    le_bare = series(args.bare_dir, 'latent_heat_flux')
    print(f"bare: LE at step 1 {le_bare[1]:.4f} W/m^2 (mid, linear factor, {le_mid[1]:.4f})")
    for n in range(1, args.steps + 1):
        if not 0.0 < le_bare[n] < le_mid[n]:
            failures.append(f"bare step {n}: LE {le_bare[n]:.4f} W/m^2 is not between 0 and the "
                            f"linear factor's {le_mid[n]:.4f}: the soil resistance is not used")
            break
    check_budget(args.bare_dir, 'bare', args.mid_q)

    # 7
    worst = 0.0
    for variable in ('q_surf', 'latent_heat_flux'):
        a = series(args.bare_dir, variable)
        b = series(args.bare_f0_dir, variable)
        for n in range(1, args.steps + 1):
            scale = max(abs(a[n]), 1.0e-30)
            worst = max(worst, abs(a[n] - b[n]) / scale)
    print(f"bare_f0 vs bare: worst relative difference of q_surf and LE {worst:.3e}")
    if worst > args.continuity_rtol:
        failures.append(f"bare_f0 differs from bare by {worst:.3e} (relative): a vegetation type "
                        f"at fraction 0 must give bare soil's beta")

    # 9
    le_bare_dry = series(args.bare_dry_dir, 'latent_heat_flux')
    worst_dry = max(le_bare_dry[1:])
    print(f"bare_dry: largest LE over the steps {worst_dry:.4f} W/m^2")
    if worst_dry > 0.0:
        failures.append(f"bare_dry: LE reaches {worst_dry:.4f} W/m^2; bare soil at the wilting "
                        f"point must not evaporate (its pore air is nearly dry)")
    check_budget(args.bare_dry_dir, 'bare_dry', args.wilt, expect_drying=False)

    # 8
    def recorded_z0(directory):
        path = os.path.join(directory, 'plt2d00000', 'job_info')
        values = []
        with open(path) as handle:
            for line in handle:
                fields = line.split('=')
                if len(fields) == 2 and fields[0].strip() == 'erf.most.z0':
                    values.append(float(fields[1].split()[0]))
        if len(values) != 1:
            raise CheckError(f"{path}: {len(values)} erf.most.z0 entries, expected one")
        return values[0]

    for label, directory, expected in (('tables', args.tables_dir, args.tables_z0),
                                       ('veg', args.veg_dir, args.given_z0)):
        z0 = recorded_z0(directory)
        print(f"{label}: job_info erf.most.z0 = {z0} (expected {expected})")
        if abs(z0 - expected) > 1.0e-6 * expected:
            failures.append(f"{label}: job_info records erf.most.z0 = {z0}, but the run used "
                            f"{expected}")

    if failures:
        for message in failures:
            print(f"FAIL: {message}")
        return 1
    print("PASS: the surface evaporates at beta times the potential rate and the soil "
          "loses the water the air gains")
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except CheckError as error:
        print(f"FAIL: {error}")
        sys.exit(1)
