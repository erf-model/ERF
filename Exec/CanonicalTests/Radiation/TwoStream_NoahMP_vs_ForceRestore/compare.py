#!/usr/bin/env python3
"""Compare the Noah-MP and force-restore runs of this case, level by level.

Reads every plt2d file of both runs and the radiation CSV of each, and for each level and
output time takes the level's mean of:

  skin temperature   Noah-MP t_sfc               force-restore seb_t_sfc
  H, LE              sensible_heat_flux, latent_heat_flux (what the surface layer applies)
  absorbed SW        Noah-MP sav + sag           force-restore SW_surface (radiation CSV)
  net LW, down       Noah-MP -fira               force-restore -LW_net_surface (radiation CSV)
  PBL height         pblh

Writes comparison.csv (one row per level and time) and comparison.png (one panel per
quantity, both runs, both levels), and prints the differences at the last output and at
the time of peak H.

Usage: compare.py --noahmp DIR --force-restore DIR [--out DIR]
"""

import argparse
import csv
import glob
import os
import sys

import numpy as np

import yt

yt.set_log_level(40)


def level_means(plotfile, names):
    """{level: {name: mean over that level's cells}} and the plotfile's time."""
    ds = yt.load(plotfile)
    out = {}
    for lev in range(ds.index.max_level + 1):
        grids = [g for g in ds.index.grids if g.Level == lev]
        out[lev] = {}
        for name in names:
            values = np.concatenate([np.asarray(g['boxlib', name]).ravel() for g in grids])
            out[lev][name] = float(values.mean())
    return float(ds.current_time), out


def radiation_csv(path):
    """{level: [(time, row)]} of the pre_dycore rows, in time order."""
    rows = {}
    with open(path) as handle:
        for r in csv.DictReader(handle):
            if r['call_site'] == 'pre_dycore':
                rows.setdefault(int(r['level']), []).append((float(r['time']), r))
    return rows


def row_at(rows, level, time, tol=5.0):
    """The level's pre_dycore row nearest to time (its sweep of that state), within tol
    seconds. Levels number their steps separately when they subcycle, so the CSV is matched
    by time, not by step."""
    best = min(rows.get(level, []), key=lambda tr: abs(tr[0] - time), default=None)
    return best[1] if best and abs(best[0] - time) <= tol else None


def plotfiles(directory):
    return sorted(glob.glob(os.path.join(directory, 'plt2d[0-9]*')),
                  key=lambda p: int(os.path.basename(p)[5:]))


def series(directory, model):
    """List of {time, step, level, skin, H, LE, SW_abs, LW_net, pblh}."""
    csv_path = glob.glob(os.path.join(directory, 'radiation_diag_*.csv'))
    if not csv_path:
        sys.exit(f"no radiation CSV in {directory}")
    rad = radiation_csv(csv_path[0])
    if model == 'noahmp':
        names = ['t_sfc', 'sav', 'sag', 'fira', 'sensible_heat_flux', 'latent_heat_flux', 'pblh']
    else:
        names = ['seb_t_sfc', 'sensible_heat_flux', 'latent_heat_flux', 'pblh']
    out = []
    for pf in plotfiles(directory):
        step = int(os.path.basename(pf)[5:])
        if step == 0:
            continue  # nothing has run yet
        time, means = level_means(pf, names)
        for lev, m in means.items():
            row = {'time_h': time / 3600.0, 'step': step, 'level': lev,
                   'H': m['sensible_heat_flux'], 'LE': m['latent_heat_flux'], 'pblh': m['pblh']}
            if model == 'noahmp':
                row.update(skin=m['t_sfc'], SW_abs=m['sav'] + m['sag'], LW_net=-m['fira'])
            else:
                r = row_at(rad, lev, time)
                # LW_net_surface is up minus down; report it positive down, as -fira is.
                row.update(skin=m['seb_t_sfc'],
                           SW_abs=float(r['SW_surface']) if r else float('nan'),
                           LW_net=-float(r['LW_net_surface']) if r else float('nan'))
            out.append(row)
    if not out:
        sys.exit(f"no plt2d output after step 0 in {directory}")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--noahmp', required=True)
    parser.add_argument('--force-restore', required=True)
    parser.add_argument('--out', default='.')
    args = parser.parse_args()

    runs = {'noahmp': series(args.noahmp, 'noahmp'),
            'force_restore': series(args.force_restore, 'force_restore')}
    keys = ['skin', 'H', 'LE', 'SW_abs', 'LW_net', 'pblh']
    labels = {'skin': 'skin temperature [K]', 'H': 'sensible heat flux H [W/m$^2$]',
              'LE': 'latent heat flux LE [W/m$^2$]', 'SW_abs': 'absorbed shortwave [W/m$^2$]',
              'LW_net': 'net longwave, down [W/m$^2$]', 'pblh': 'PBL height [m]'}

    with open(os.path.join(args.out, 'comparison.csv'), 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['model', 'level', 'time_h', 'step'] + keys)
        for model, rows in runs.items():
            for r in rows:
                writer.writerow([model, r['level'], f"{r['time_h']:.4f}", r['step']]
                                + [f"{r[k]:.6g}" for k in keys])

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    styles = {('noahmp', 0): ('#1f77b4', '-'), ('noahmp', 1): ('#1f77b4', '--'),
              ('force_restore', 0): ('#d62728', '-'), ('force_restore', 1): ('#d62728', '--')}
    for ax, key in zip(axes.ravel(), keys):
        for model, rows in runs.items():
            for lev in (0, 1):
                sel = [r for r in rows if r['level'] == lev]
                color, ls = styles[(model, lev)]
                ax.plot([r['time_h'] for r in sel], [r[key] for r in sel], color=color, ls=ls,
                        label=f"{'Noah-MP' if model == 'noahmp' else 'force-restore'}, level {lev}")
        ax.set_title(labels[key], fontsize=10)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel('hours after 15:00 UTC')
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'comparison.png'), dpi=120)

    # The two runs take different time steps (each follows its own CFL), so their outputs
    # carry different step numbers; pair them by output time (every 10 minutes).
    def slot(r):
        return round(r['time_h'] * 6.0)

    def at(model, lev, s):
        return next((r for r in runs[model] if r['level'] == lev and slot(r) == s), None)

    last = max(slot(r) for r in runs['noahmp'])
    peak = slot(max((r for r in runs['noahmp'] if r['level'] == 0), key=lambda r: r['H']))
    for title, s in (('last output', last), ('peak Noah-MP H', peak)):
        print(f"{title} ({s / 6.0:.2f} h after 15:00 UTC):")
        for lev in (0, 1):
            n, f = at('noahmp', lev, s), at('force_restore', lev, s)
            if n is None or f is None:
                print(f"  level {lev}: no output at that time in both runs")
                continue
            print(f"  level {lev}: " + ', '.join(
                f"{k} {n[k]:.4g} vs {f[k]:.4g}" for k in keys))
    return 0


if __name__ == '__main__':
    sys.exit(main())
