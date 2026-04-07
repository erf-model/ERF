import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from erftools.postprocessing import Plotfile

plt.rcParams['axes.labelsize'] = 'x-large'
plt.rcParams['xtick.labelsize'] = 'large'
plt.rcParams['ytick.labelsize'] = 'large'
plt.rcParams['font.family'] = 'serif'

def animate_pltfiles(simdir):
    animdir = f'{simdir}/anim'
    os.makedirs(animdir, exist_ok=True)

    pltfiles = glob.glob(f'{simdir}/plt*')
    pltfiles = [dpath for dpath in pltfiles if os.path.isdir(dpath) and not '.old.' in dpath]
    outsteps = [int(os.path.split(dpath)[1][3:]) for dpath in pltfiles]
    pltorder = np.argsort(outsteps)
    pltfiles = [pltfiles[idx] for idx in pltorder]

    dpath = pltfiles[0]
    pf = Plotfile(dpath)

    yslice = pf.prob_extent[1] / 2

    # setup plotting grid
    nx,ny,nz = pf.dims
    x1 = np.linspace(0, pf.prob_extent[0], nx+1)
    if pf.stretched_dz:
        z1 = pf.terrain_z_levels
    else:
        z1 = np.linspace(0, pf.prob_extent[2], nz+1)
    xx,zz = np.meshgrid(x1, z1, indexing='ij')

    # plot snapshots
    fig,ax = plt.subplots(figsize=(13,4))
    for itime, dpath in enumerate(pltfiles):
        if itime > 0:
            pf = Plotfile(dpath)
        ds = pf.slice(1, yslice, fields='x_velocity')
        fld = ds['x_velocity'].squeeze()
        print(f'\r({itime+1}/{len(pltfiles)}) t={pf.time:g}',end='')

        if itime == 0:
            cmsh = ax.pcolormesh(xx, zz, fld, vmin=0)
            cbar = fig.colorbar(cmsh, ax=ax, label='$u^+$')
            ax.set_xlabel('$x/\\delta$')
            ax.set_ylabel('$y/\\delta$')
            ax.axis('scaled')
        else:
            cmsh.set_array(fld.values.ravel())
            #cmsh.set_clim((np.min(values), np.max(values)))

        ax.set_title(f'$t = {pf.time:.4f}$')

        fig.savefig(f'{animdir}/channel_DNS_u.{itime:04d}.png')

#==============================================================================
if __name__ == '__main__':
    simdir = sys.argv[1]
    animate_pltfiles(simdir)
