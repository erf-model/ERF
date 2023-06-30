#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yt
from couette_helpers import read_input

#outputs = [100, 200, 1300, 2600, 12800]
outputs = [10,20,30,40,50]
Utop = 2.0

try:
    results_dir = sys.argv[1]
except IndexError:
    results_dir = '.'
pltfiles = [os.path.join(results_dir,f'plt{idx:05d}') for idx in outputs]
assert all([os.path.isdir(pltfile) for pltfile in pltfiles])

pp = read_input(os.path.join(results_dir,'inputs_couette'))
dt = pp['erf.fixed_dt']
print('dt=',dt)
nu = pp['erf.dynamicViscosity'] / pp['prob.rho_0']
print('nu=',nu)
zmin = float(pp['geometry.prob_lo'][2])
zmax = float(pp['geometry.prob_hi'][2])
nz = int(pp['amr.n_cell'][2])
h = zmax - zmin
print('h=',h)

zface = np.linspace(zmin,zmax,nz+1)
z = (zface[1:] + zface[:-1]) / 2

nondimtimes = nu * np.array(outputs)*dt / h**2

#
# plot evolution of profiles
#
fig,ax = plt.subplots()
for ndt,pltfile in zip(nondimtimes,pltfiles):
    print(ndt,pltfile)
    pf = yt.load(os.path.join(results_dir, pltfile))
    level0 = pf.covering_grid(level=0, left_edge=pf.domain_left_edge, dims=pf.domain_dimensions)
    u = level0['x_velocity'][0,:,0]
    ax.plot(u/Utop, z/h, label=ndt)
ax.set_xlim((0,1))
ax.set_ylim((0,1))
lgd = ax.legend(title=r'$\nu t/h^2$',title_fontsize='large')
ax.set_ylabel(r'$z/h$',fontsize='x-large')
ax.set_xlabel(r'$u/U$',fontsize='x-large')
ax.grid()
fig.savefig('Couette_profiles.png',bbox_inches='tight')

#
# check error at last step
#
err = u.value - Utop*z/h
print('MAE',np.mean(np.abs(err)))
print('RMSE', np.sqrt(np.mean(err*err)))
