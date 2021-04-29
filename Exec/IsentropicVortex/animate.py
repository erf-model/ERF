#!/usr/bin/env python
import glob
import numpy as np
import matplotlib.pyplot as plt
import yt
from isentropic_vortex_helpers import read_input

cmaps = {
    'density':'viridis',
    'temp':'magma',
    'x_velocity':'RdBu_r',
    'y_velocity':'RdBu_r',
}

prob_parm = read_input('inputs_ex',prefix='prob.')
print('Problem parameters from inputs_ex:\n',prob_parm)

ainf = np.sqrt(prob_parm['gamma']*287.0*prob_parm['T_inf'])
uinf = prob_parm['M_inf'] * ainf * np.cos(prob_parm['alpha'])
vinf = prob_parm['M_inf'] * ainf * np.sin(prob_parm['alpha'])
print('speed of sound:',ainf)
print('freestream:',uinf,vinf)

pltfiles = glob.glob('plt?????')
pltfiles.sort()
print('Solutions found:',pltfiles)

xx, yy = None, None
kidx = 0
fields = ['density','temp','x_velocity','y_velocity']
minmax = {field:None for field in fields} # plot ranges
for dpath in pltfiles:
    print(f'Processing {dpath}')
    ds = yt.load(dpath)
    if xx is None:
        # get grid points once--no AMR
        pts1d = [
        np.linspace(float(x0),float(x1),N+1)
        for x0,x1,N in zip(ds.domain_left_edge,
                           ds.domain_right_edge,
                           ds.domain_dimensions)
        ]
        xx,yy,_ = np.meshgrid(*pts1d,indexing='ij')
        xx = xx[:,:,kidx]
        yy = yy[:,:,kidx]

    # new figure (less efficient for now)
    fig,ax = plt.subplots(nrows=2,ncols=2,
                          sharex=True,sharey=True,
                          figsize=(8,8))
    tsim = float(ds.current_time)
    fig.suptitle(f't = {tsim:g} s')

    # get data (only level 0, no AMR)
    level0 = ds.covering_grid(level=0,
                              left_edge=ds.domain_left_edge,
                              dims=ds.domain_dimensions)

    # plot snapshots
    for axi,field in zip(ax.flatten(),fields):
        fielddata = level0[field][:,:,kidx].value
        if field == 'x_velocity':
            fielddata -= uinf
        if field == 'y_velocity':
            fielddata -= vinf
        if minmax[field] is None:
            minmax[field] = dict(vmin=np.min(fielddata), vmax=np.max(fielddata))
#        plotparams = {'cmap':cmaps[field]}
#        if field.endswith('velocity'):
#            # make symmetric about 0
#            plotparams['vmin'] = min(np.min(fielddata), -np.max(fielddata))
#            plotparams['vmax'] = max(np.max(fielddata), -np.min(fielddata))
#        cm = axi.pcolormesh(xx,yy,fielddata,**plotparams)
        cm = axi.pcolormesh(xx,yy,fielddata,cmap=cmaps[field],**minmax[field])
        cb = fig.colorbar(cm,ax=axi)
        cb.set_label(field,fontsize='x-large')
        axi.axis('scaled')

    fig.savefig(f'{dpath}.png',dpi=150)
    plt.close()

