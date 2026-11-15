#!/usr/bin/env python3
"""Check the force-restore snow layer of the two-stream surface energy balance.

The runs are the TwoStream_SEBSnow deck with the 2D plotfile written every step
(Tests/RunTwoStreamSEBSnow.cmake). The closed forms are those of Douville, Royer and Mahfouf
(1995, Climate Dynamics 12, 21-35; "DRM95"), Sec. 2.1. Every check is made per cell.

melt (a 50 kg/m^2 snowpack under a warm sky; it does not melt away), every step n >= 1:

1. the surface is at the melting point, T_s(n) = 273.16 K (DRM95 Eq. 12 and Sec. 2.1);
2. the snow budget closes (DRM95 Eq. 6): W_n(n) - W_n(0) = snowfall - sublimation - melt, with
   the running totals from the plotfile, and the snow has melted in every step;
3. the step's sublimation is the covered fraction of the evaporation (DRM95 Eq. 7):
     sublimation(n) - sublimation(n-1) = dt max(LE(n), 0) / L_v * W_n(n-1) / (W_n(n-1) + W_crn);
4. the absorbed shortwave is that of the albedo with snow (DRM95 Eq. 8 with veg = 0),
     SW_net(n) = SW_net_bare(n) (1 - alpha_t) / (1 - alpha_g),
     alpha_t = alpha_g + (alpha_n - alpha_g) W_n(n-1) / (W_n(n-1) + W_crn),
   against the bare leg (the same deck without snow), whose downwelling flux is the same.

full_melt (a 1 kg/m^2 snowpack with the snow albedo set to the ground's, so the first step's
balance is the bare leg's), at step 1:

5. all the snow is gone, and melt + sublimation = 1 kg/m^2;
6. the surface is colder than the bare leg's by the heat of the melt (DRM95 Eq. 12),
     T_s_bare(1) - T_s(1) = L_f melt / C_s,
   and still above the melting point;
7. the soil has the snow's water: q_s(1) - q_s_bare(1) = (melt + sublimation) / (rho_w d_s).

two_level (the melt leg with a second level made half-way through): checks 1 and 2 on every
cell of the finest data at the last step, the snow and the melt total there equal the melt
leg's (the new level is filled from the coarse one), and the last plotfile has level 1.
"""

import argparse
import glob
import os
import subprocess
import sys

L_V = 2.5e6       # latent heat of vaporization [J/kg] (ERF_Constants.H)
L_F = 3.337e5     # latent heat of fusion [J/kg] (lat_ice, ERF_MicrophysicsConstants.H)
RHO_W = 1000.0    # density of water [kg/m^3]
T_00 = 273.16     # melting point of DRM95 Eq. (12) [K]


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


def on_cells(points, xs):
    """The value of POINTS at each x of XS: the cell centre nearest x (the coarse cell that
    contains x where the plotfile is coarser than XS)."""
    return [min(points, key=lambda p: abs(p[0] - x))[1] for x in xs]


def plotfiles(run_dir, steps):
    files = sorted(f for f in glob.glob(os.path.join(run_dir, 'plt2d*')) if os.path.isdir(f))
    if len(files) != steps + 1:
        raise CheckError(f"{run_dir}: expected {steps + 1} 2D plotfiles (steps 0..{steps}), "
                         f"found {len(files)}")
    return files


def close(a, b, rtol, atol, what):
    if not abs(a - b) <= atol + rtol * max(abs(a), abs(b)):
        raise CheckError(f"{what}: {a!r} vs {b!r} (difference {a - b:.3e})")


def finest_level(plotfile):
    with open(os.path.join(plotfile, 'Header')) as handle:
        lines = handle.read().split('\n')
    nvars = int(lines[1])
    return int(lines[2 + nvars + 2])  # after the names: dimension, time, finest level


class Run:
    """The plotfiles of one leg, read on the cells of its last plotfile."""

    def __init__(self, args, run_dir, steps):
        self.args = args
        self.files = plotfiles(run_dir, steps)
        self.tmp = os.path.join(run_dir, 'fextract.tmp')
        self.xs = [x for x, _ in cell_points(args.fextract, self.files[-1], 'seb_t_sfc', self.tmp)]

    def get(self, n, variable):
        return on_cells(cell_points(self.args.fextract, self.files[n], variable, self.tmp), self.xs)


def fraction(w, w_crn):
    return w / (w + w_crn) if w > 0.0 else 0.0


def check_budget(run, n, w0, label):
    swe, fall, sub, melt = (run.get(n, v) for v in
                            ('seb_snow_swe', 'seb_snowfall_total',
                             'seb_snow_sublimation_total', 'seb_snow_melt_total'))
    for c in range(len(swe)):
        close(swe[c] - w0[c], fall[c] - sub[c] - melt[c], 1e-12, 1e-11,
              f"{label}: snow budget at step {n}, cell {c} (DRM95 Eq. 6)")
    return swe, fall, sub, melt


def check_melt(args):
    run = Run(args, args.melt_dir, args.steps)
    bare = Run(args, args.bare_dir, args.steps)
    w0 = run.get(0, 'seb_snow_swe')
    for c, w in enumerate(w0):
        close(w, args.swe0, 1e-14, 0.0, f"melt: seb_snow_swe at step 0, cell {c}")
    w_prev, sub_prev, melt_prev = w0, [0.0] * len(w0), [0.0] * len(w0)
    for n in range(1, args.steps + 1):
        swe, _, sub, melt = check_budget(run, n, w0, 'melt')
        ts, lh, sw = (run.get(n, v) for v in ('seb_t_sfc', 'seb_lh', 'seb_sw_net'))
        sw_bare = bare.get(n, 'seb_sw_net')
        for c in range(len(swe)):
            if not swe[c] > 0.0:
                raise CheckError(f"melt: the snow melted away at step {n}, cell {c}; the "
                                 "checks below assume it stays")
            close(ts[c], T_00, 1e-13, 0.0, f"melt: T_s at step {n}, cell {c} (DRM95 Eq. 12)")
            if not melt[c] > melt_prev[c]:
                raise CheckError(f"melt: no melt in step {n}, cell {c} "
                                 f"({melt_prev[c]} -> {melt[c]})")
            p_ng = fraction(w_prev[c], args.w_crn)
            if not (lh[c] > 0.0 and sub[c] > sub_prev[c]):
                raise CheckError(f"melt: no sublimation in step {n}, cell {c} (LE = {lh[c]}); "
                                 "the check below needs an upward flux")
            close(sub[c] - sub_prev[c], args.dt * max(lh[c], 0.0) / L_V * p_ng, 1e-7, 1e-13,
                  f"melt: sublimation of step {n}, cell {c} (DRM95 Eq. 7)")
            alpha_t = args.albedo_ground + (args.albedo_snow - args.albedo_ground) * p_ng
            if not sw_bare[c] > 1.0:
                raise CheckError(f"melt: the bare leg absorbs no shortwave at step {n}")
            close(sw[c], sw_bare[c] * (1.0 - alpha_t) / (1.0 - args.albedo_ground), 1e-9, 0.0,
                  f"melt: absorbed shortwave at step {n}, cell {c} (DRM95 Eq. 8)")
        w_prev, sub_prev, melt_prev = swe, sub, melt
    print(f"PASS melt: {args.steps} steps at {T_00} K, snow {args.swe0} -> {min(w_prev):.6f} "
          f"kg/m^2, melt {max(melt_prev):.6e}, sublimation {max(sub_prev):.6e} kg/m^2; budget, "
          "sublimation and albedo per cell")


def check_full_melt(args):
    run = Run(args, args.full_melt_dir, args.steps)
    bare = Run(args, args.bare_dir, args.steps)
    w0 = run.get(0, 'seb_snow_swe')
    swe, _, sub, melt = check_budget(run, 1, w0, 'full_melt')
    ts, qs = run.get(1, 'seb_t_sfc'), run.get(1, 'seb_q_sfc')
    ts_bare, qs_bare = bare.get(1, 'seb_t_sfc'), bare.get(1, 'seb_q_sfc')
    for c in range(len(swe)):
        close(w0[c], args.full_melt_swe0, 1e-14, 0.0, f"full_melt: seb_snow_swe at step 0, cell {c}")
        if swe[c] != 0.0:
            raise CheckError(f"full_melt: {swe[c]} kg/m^2 of snow left after step 1, cell {c}")
        close(melt[c] + sub[c], args.full_melt_swe0, 1e-12, 0.0,
              f"full_melt: melt + sublimation, cell {c}")
        close(ts_bare[c] - ts[c], L_F * melt[c] / args.c_s, 1e-9, 0.0,
              f"full_melt: cooling by the melt, cell {c} (DRM95 Eq. 12)")
        if not ts[c] > T_00:
            raise CheckError(f"full_melt: T_s {ts[c]} not above the melting point with no snow "
                             f"left, cell {c}")
        close(qs[c] - qs_bare[c], (melt[c] + sub[c]) / (RHO_W * args.d_s), 1e-9, 0.0,
              f"full_melt: snow water in the soil, cell {c}")
    print(f"PASS full_melt: {args.full_melt_swe0} kg/m^2 melted in one step, the surface "
          f"{max(ts_bare) - max(ts):.6f} K colder than bare ground (L_f melt / C_s), the soil "
          "holds the water")


def check_two_level(args):
    run = Run(args, args.two_level_dir, args.steps)
    last = args.steps
    if finest_level(run.files[last]) != 1:
        raise CheckError(f"two_level: {run.files[last]} has no level 1")
    w0 = run.get(0, 'seb_snow_swe')
    swe, _, _, melt = check_budget(run, last, w0, 'two_level')
    ts = run.get(last, 'seb_t_sfc')
    # The column is the melt leg's and the flow is uniform in x and y, so a level made
    # half-way must carry the melt leg's snow and totals (filled from the coarse level), not
    # the snow of the start.
    one = Run(args, args.melt_dir, args.steps)
    swe_one = one.get(last, 'seb_snow_swe')[0]
    melt_one = one.get(last, 'seb_snow_melt_total')[0]
    for c in range(len(swe)):
        if not (swe[c] > 0.0 and melt[c] > 0.0):
            raise CheckError(f"two_level: snow {swe[c]}, melt {melt[c]} at cell {c}")
        close(ts[c], T_00, 1e-13, 0.0, f"two_level: T_s at step {last}, cell {c}")
        close(swe[c], swe_one, 1e-6, 0.0,
              f"two_level: seb_snow_swe at step {last}, cell {c}, against the melt leg")
        close(melt[c], melt_one, 1e-6, 0.0,
              f"two_level: seb_snow_melt_total at step {last}, cell {c}, against the melt leg")
    print(f"PASS two_level: {len(swe)} cells on the finest data at {T_00} K, the budget "
          "closes per cell")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fextract', required=True)
    parser.add_argument('--melt-dir', required=True)
    parser.add_argument('--bare-dir', required=True)
    parser.add_argument('--full-melt-dir', required=True)
    parser.add_argument('--two-level-dir', required=True)
    parser.add_argument('--steps', type=int, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--swe0', type=float, required=True, help='initial SWE of melt [kg/m^2]')
    parser.add_argument('--full-melt-swe0', type=float, required=True)
    parser.add_argument('--w-crn', type=float, required=True, help='seb_snow_critical_swe_kg_m2')
    parser.add_argument('--albedo-snow', type=float, required=True, help='seb_snow_albedo')
    parser.add_argument('--albedo-ground', type=float, required=True, help='surface_albedo_sw')
    parser.add_argument('--c-s', type=float, required=True, help='seb_surface_heat_capacity')
    parser.add_argument('--d-s', type=float, required=True, help='seb_moisture_layer_depth_m')
    args = parser.parse_args()
    checks = [check_melt, check_full_melt, check_two_level]
    failed = False
    for check in checks:
        try:
            check(args)
        except CheckError as err:
            print(f"FAIL {err}")
            failed = True
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
