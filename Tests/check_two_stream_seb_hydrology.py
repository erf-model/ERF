#!/usr/bin/env python3
"""Check the force-restore hydrology of the two-stream surface energy balance.

The runs are the TwoStream_SEBSurfaceLayerFluxes deck (Kessler, surface-layer fluxes) with the
2D plotfile written every step. Given the plotfiles of two runs this asserts the following.

deep (seb_deep_prognostic_enable and seb_drainage_enable, no rain), every step n >= 1, on the
column mean (all columns stay above field capacity, so the relations are linear in the state):

1. the column water follows MN96 Eq. (4) with E = LE / L_v, the drainage exact over the step:
     q_deep(n) = w_fc + (q_deep(n-1) - dt LE(n) / (L_v rho_w d_2) - w_fc) exp(-C_3 dt / tau_1d),
   with LE(n) the latent heat flux the balance used (seb_lh) and tau_1d = 86400 s, the period
   of MN96 Table 1 (not seb_restore_timescale_s), and it has drained;
2. the drainage total grows by the water that left: rho_w d_2 times the drop above;
3. the deep temperature relaxes exactly toward the updated skin temperature (NP89 Eq. 6):
     T_deep(n) = T_s(n) + (T_deep(n-1) - T_s(n)) exp(-dt / tau), tau = seb_restore_timescale_s.

rain (seb_precipitation_source = microphysics with the deep reservoir and drainage, a moist
sounding so that Kessler rains, and a saturated soil, q_s = q_deep = porosity), at the last step:

4. rain reached the ground, and the rain the soil took in equals the change of the scheme's
   own accumulation, seb_rain_total = precip_total_accum(n) - precip_total_accum(0), per cell;
5. the column budget closes,
     rho_w d_2 (q_deep(n) - q_deep(0)) = rain - sum(E dt) - drainage - runoff + deficit,
   with the totals from the plotfile and E from seb_lh at every step;
6. water ran off (rain on a saturated soil, once it exceeds the evaporation), and the surface
   water content does not exceed the porosity.

two_level (the rain leg with a second level made half-way through, regridded every step, with
no subcycling): checks 4-6 on every column of the finest data at the last step, and the last
plotfile has level 1. Earlier plotfiles are read on the last one's cells (a fine cell takes the
coarse cell that contains it before level 1 exists).

two_level_subcycled (the two_level run with level 1 taking two steps per coarse step, so the
coarse totals and accumulation are averaged down from level 1 between coarse steps): checks 4
and 6 on the finest data. Check 5 is not made: the plotfile holds only the last substep's latent
heat flux, so the evaporation of level 1 cannot be summed from it.
"""

import argparse
import glob
import math
import os
import subprocess
import sys

L_V = 2.5e6       # latent heat of vaporization [J/kg] (ERF_Constants.H)
RHO_W = 1000.0    # density of water [kg/m^3] (rhor, ERF_MicrophysicsConstants.H)
TAU_C3 = 86400.0  # restore period of MN96 Table 1 [s] (force_restore_c3_tau)


class CheckError(Exception):
    """A condition that fails the check (reported as FAIL, not as a traceback)."""


def cell_points(fextract, plotfile, variable, out_path):
    """(x, value) of every cell on the x slice, from the finest data covering it."""
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
                values.append((float(fields[0]), float(fields[1])))
    if not values:
        raise CheckError(f"{plotfile}: no {variable} cells on the slice")
    return values


def cell_values(fextract, plotfile, variable, out_path):
    return [v for _, v in cell_points(fextract, plotfile, variable, out_path)]


def on_cells(points, xs):
    """The value of POINTS at each x of XS: the cell centre nearest x. Where a plotfile is
    coarser than XS this is the coarse cell that contains x (centres of a coarser grid are at
    least half a coarse cell from a fine centre in a neighbouring coarse cell)."""
    return [min(points, key=lambda p: abs(p[0] - x))[1] for x in xs]


def mean(values):
    return sum(values) / len(values)


def plotfiles(run_dir, steps):
    files = sorted(glob.glob(os.path.join(run_dir, 'plt2d*')))
    files = [f for f in files if os.path.isdir(f)]
    if len(files) != steps + 1:
        raise CheckError(f"{run_dir}: expected {steps + 1} 2D plotfiles (steps 0..{steps}), "
                         f"found {len(files)}")
    return files


def close(a, b, rtol, atol, what):
    if not abs(a - b) <= atol + rtol * max(abs(a), abs(b)):
        raise CheckError(f"{what}: {a!r} vs {b!r} (difference {a - b:.3e})")


def check_deep(args):
    files = plotfiles(args.deep_dir, args.deep_steps)
    tmp = os.path.join(args.deep_dir, 'fextract.tmp')
    get = lambda f, v: mean(cell_values(args.fextract, f, v, tmp))
    k_drain = math.exp(-args.c3_m / args.d2 * args.dt / TAU_C3)
    k_temp = math.exp(-args.dt / args.tau)
    qd_prev = get(files[0], 'seb_q_deep')
    td_prev = get(files[0], 'seb_t_deep')
    dr_prev = get(files[0], 'seb_drainage_total')
    close(qd_prev, args.q_deep0, 1e-14, 0.0, "deep: q_deep at step 0")
    for n in range(1, len(files)):
        f = files[n]
        lh, qd, td, ts, dr = (get(f, v) for v in
                              ('seb_lh', 'seb_q_deep', 'seb_t_deep', 'seb_t_sfc',
                               'seb_drainage_total'))
        pre = qd_prev - args.dt * lh / (L_V * RHO_W * args.d2)
        if not pre > args.fc:
            raise CheckError(f"deep: the column fell to field capacity at step {n}; "
                             "the closed form below assumes it stays above")
        expect = args.fc + (pre - args.fc) * k_drain
        close(qd, expect, 1e-12, 0.0, f"deep: q_deep at step {n} (MN96 Eq. 4)")
        close(dr - dr_prev, RHO_W * args.d2 * (pre - qd), 1e-8, 1e-11,
              f"deep: drainage at step {n}")
        close(td, ts + (td_prev - ts) * k_temp, 1e-13, 0.0,
              f"deep: T_deep at step {n} (NP89 Eq. 6)")
        qd_prev, td_prev, dr_prev = qd, td, dr
    if not dr_prev > 0.0:
        raise CheckError(f"deep: the column did not drain (drainage total {dr_prev})")
    print(f"PASS deep: {len(files) - 1} steps, q_deep {args.q_deep0} -> {qd_prev:.10f}, "
          f"drainage {dr_prev:.6e} kg/m^2")


def finest_level(plotfile):
    with open(os.path.join(plotfile, 'Header')) as handle:
        lines = handle.read().split('\n')
    nvars = int(lines[1])
    return int(lines[2 + nvars + 2])  # after the names: dimension, time, finest level


def check_rain(args, run_dir=None, label='rain', budget=True):
    run_dir = run_dir or args.rain_dir
    files = plotfiles(run_dir, args.rain_steps)
    tmp = os.path.join(run_dir, 'fextract.tmp')
    cells = lambda f, v: cell_values(args.fextract, f, v, tmp)
    last = files[-1]
    # The cells of the last plotfile; earlier plotfiles (coarser where level 1 did not yet
    # exist) are read on them.
    xs = [x for x, _ in cell_points(args.fextract, last, 'seb_q_deep', tmp)]
    on_last = lambda f, v: on_cells(cell_points(args.fextract, f, v, tmp), xs)
    acc0 = on_last(files[0], 'precip_total_accum')
    qd0 = on_last(files[0], 'seb_q_deep')
    evap = [0.0] * len(xs)
    for f in files[1:]:
        lh = on_last(f, 'seb_lh')
        evap = [e + args.dt * l / L_V for e, l in zip(evap, lh)]
    acc = cells(last, 'precip_total_accum')
    rain = cells(last, 'seb_rain_total')
    runoff = cells(last, 'seb_runoff_total')
    drain = cells(last, 'seb_drainage_total')
    deficit = cells(last, 'seb_deficit_total')
    qd = cells(last, 'seb_q_deep')
    qs = cells(last, 'seb_q_sfc')
    if len(rain) != len(xs):
        raise CheckError(f"{label}: {len(xs)} cells of seb_q_deep but {len(rain)} of "
                         "seb_rain_total in the last plotfile")
    if not max(rain) > 0.0:
        raise CheckError(f"{label}: no rain reached the ground in {args.rain_steps} steps")
    for c in range(len(rain)):
        close(rain[c], acc[c] - acc0[c], 1e-10, 1e-14, f"{label}: seb_rain_total, cell {c}")
        if budget:
            lhs = RHO_W * args.d2 * (qd[c] - qd0[c])
            rhs = rain[c] - evap[c] - drain[c] - runoff[c] + deficit[c]
            close(lhs, rhs, 1e-9, 1e-12, f"{label}: column budget, cell {c}")
        if qs[c] > args.porosity * (1.0 + 1e-14):
            raise CheckError(f"{label}: q_s {qs[c]} above the porosity {args.porosity}, cell {c}")
    if not max(runoff) > 0.0:
        raise CheckError(f"{label}: no runoff, although rain fell on a saturated soil")
    extra = ''
    if label.startswith('two_level'):
        if finest_level(last) != 1:
            raise CheckError(f"{label}: {last} has no level 1")
        extra = f', {len(rain)} cells on the finest data'
    print(f"PASS {label}: {args.rain_steps} steps, mean rain {mean(rain):.6e}, runoff "
          f"{mean(runoff):.6e}, drainage {mean(drain):.6e} kg/m^2"
          + ("; the budget closes per cell" if budget else "; rain matches the accumulation per cell")
          + extra)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fextract', required=True)
    parser.add_argument('--deep-dir', required=True)
    parser.add_argument('--deep-steps', type=int, required=True)
    parser.add_argument('--rain-dir', required=True)
    parser.add_argument('--rain-steps', type=int, required=True)
    parser.add_argument('--two-level-dir', required=True)
    parser.add_argument('--subcycled-dir', required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--tau', type=float, required=True, help='seb_restore_timescale_s')
    parser.add_argument('--d2', type=float, required=True, help='seb_deep_layer_depth_m')
    parser.add_argument('--c3-m', type=float, required=True, help='seb_drainage_c3_m')
    parser.add_argument('--fc', type=float, required=True, help='seb_soil_moisture_fc')
    parser.add_argument('--porosity', type=float, required=True, help='seb_soil_porosity')
    parser.add_argument('--q-deep0', type=float, required=True, help='seb_q_deep_default (deep)')
    args = parser.parse_args()
    failed = False
    checks = (check_deep, check_rain,
              lambda a: check_rain(a, a.two_level_dir, 'two_level'),
              lambda a: check_rain(a, a.subcycled_dir, 'two_level_subcycled', budget=False))
    for check in checks:
        try:
            check(args)
        except CheckError as err:
            print(f"FAIL {err}")
            failed = True
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
