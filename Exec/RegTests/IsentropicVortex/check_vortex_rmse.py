#!/usr/bin/env python
import sys
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yt
from isentropic_vortex_helpers import read_input, IsentropicVortex

try:
    results_dir = sys.argv[1]
except IndexError:
    results_dir = '.'

prob_parm = read_input(os.path.join(results_dir,'inputs'),prefix='prob.')

pltfiles = glob.glob(os.path.join(results_dir,'plt??????'))
pltfiles.sort()
print('Solutions found:',pltfiles)

ds = yt.load(os.path.join(results_dir, pltfiles[0]))

# generate reference solution
refsoln = IsentropicVortex(**prob_parm)
ent0 = refsoln.p_inf / refsoln.rho_inf**refsoln.gamma

# calculate rms errors
it = []
t = []
rmse_ent = [] # p/rho**gamma
rmse_rho = []
rmse_u = []
rmse_v = []
rmse_T = []
mae_ent = []
mae_rho = []
mae_u = []
mae_v = []
mae_T = []
for idx,pltfile in enumerate(pltfiles):
    if idx > 0:
        ds = yt.load(pltfile)
    itstr = os.path.split(pltfile)[1][len('plt'):]
    it.append(int(itstr))
    t.append(ds.current_time.value)
    level0 = ds.covering_grid(level=0, left_edge=[0,0.0,0.0], dims=ds.domain_dimensions)

    xx,yy,rho0,u0,v0,p0,T0 = refsoln.evaluate_at_time(
        ds.current_time.value,
        ds.domain_left_edge,
        ds.domain_right_edge,
        ds.domain_dimensions
    )

    ent = level0['pressure'].value / level0['density'].value**refsoln.gamma
    err_ent = (ent - ent0).flatten()

    err_rho = np.array([])
    err_u = np.array([])
    err_v = np.array([])
    err_T = np.array([])
    # loop over z slices
    for k in range(ds.domain_dimensions[2]):
        err_rho = np.concatenate([err_rho, (level0['density'][:,:,k].value - rho0).flatten()])
        err_u = np.concatenate([err_u, (level0['x_velocity'][:,:,k].value - u0).flatten()])
        err_v = np.concatenate([err_v, (level0['y_velocity'][:,:,k].value - v0).flatten()])
        err_T = np.concatenate([err_T, (level0['temp'][:,:,k].value - T0).flatten()])
    # RMSE calc
    rmse_ent.append(np.sqrt(np.mean(err_ent**2)))
    rmse_rho.append(np.sqrt(np.mean(err_rho**2)))
    rmse_u.append(np.sqrt(np.mean(err_u**2)))
    rmse_v.append(np.sqrt(np.mean(err_v**2)))
    rmse_T.append(np.sqrt(np.mean(err_T**2)))
    # MAE calc
    mae_ent.append(np.mean(np.abs(err_ent)))
    mae_rho.append(np.mean(np.abs(err_rho)))
    mae_u.append(np.mean(np.abs(err_u)))
    mae_v.append(np.mean(np.abs(err_v)))
    mae_T.append(np.mean(np.abs(err_T)))

# normalize
#rmse_rho = np.array(rmse_rho) / refsoln.rho_inf
#rmse_u = np.array(rmse_u) / refsoln.a_inf
#rmse_v = np.array(rmse_v) / refsoln.a_inf
#rmse_T = np.array(rmse_T) / refsoln.T_inf

# save error data
outfile = os.path.join(results_dir, 'rmse.csv')
itidx = pd.Index(it,name='iteration')
rmse = pd.DataFrame(
    {
        't':t,
        'density':rmse_rho,
        'x_velocity':rmse_u,
        'y_velocity':rmse_v,
        'temp':rmse_T,
        'p_rhogam':rmse_ent,
    },
    index=itidx
)
rmse.to_csv(outfile)
print('wrote',outfile)

outfile = os.path.join(results_dir, 'mae.csv')
mae = pd.DataFrame(
    {
        't':t,
        'density':mae_rho,
        'x_velocity':mae_u,
        'y_velocity':mae_v,
        'temp':mae_T,
        'p_rhogam':mae_ent,
    },
    index=itidx
)
mae.to_csv(outfile)
print('wrote',outfile)

# plot
fig,ax = plt.subplots()
ax.semilogy(rmse.index, rmse['density'], linestyle='none', marker='.', label='density')
ax.semilogy(rmse.index, rmse['x_velocity'], linestyle='none', marker='.', label='x-velocity')
ax.semilogy(rmse.index, rmse['y_velocity'], linestyle='none', marker='.', label='y-velocity')
ax.semilogy(rmse.index, rmse['temp'], linestyle='none', marker='.', label='temperature')
ax.set_ylim(bottom=1e-7)
ax.set_xlabel('iteration',fontsize='x-large')
ax.set_ylabel('RMSE [-]',fontsize='x-large')
ax.legend(loc='best',ncol=2)
outfile = os.path.join(results_dir, 'rmse.png')
fig.savefig(outfile,bbox_inches='tight')
print('saved',outfile)

fig,ax = plt.subplots()
ax.semilogy(mae.index, mae['density'], linestyle='none', marker='.', label='density')
ax.semilogy(mae.index, mae['x_velocity'], linestyle='none', marker='.', label='x-velocity')
ax.semilogy(mae.index, mae['y_velocity'], linestyle='none', marker='.', label='y-velocity')
ax.semilogy(mae.index, mae['temp'], linestyle='none', marker='.', label='temperature')
ax.set_ylim(bottom=1e-7)
ax.set_xlabel('iteration',fontsize='x-large')
ax.set_ylabel('MAE [-]',fontsize='x-large')
ax.legend(loc='best',ncol=2)
outfile = os.path.join(results_dir, 'mae.png')
fig.savefig(outfile,bbox_inches='tight')
print('saved',outfile)

