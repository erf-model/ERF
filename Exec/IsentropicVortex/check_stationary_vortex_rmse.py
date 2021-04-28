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

prob_parm = read_input(os.path.join(results_dir,'inputs_ex'),prefix='prob.')

pltfiles = glob.glob(os.path.join(results_dir,'plt?????'))
pltfiles.sort()
print('Solutions found:',pltfiles)

ds = yt.load(os.path.join(results_dir, pltfiles[0]))

# generate reference solution
refsoln = IsentropicVortex(**prob_parm)
xx,yy,rho0,u0,v0,p0,T0 = refsoln.evaluate(ds.domain_left_edge,
                                          ds.domain_right_edge,
                                          ds.domain_dimensions)

# calculate rms errors
it = []
t = []
rmse = []
for idx,pltfile in enumerate(pltfiles):
    if idx > 0:
        ds = yt.load(pltfile)
    itstr = os.path.split(pltfile)[1][len('plt'):]
    it.append(int(itstr))
    t.append(ds.current_time.value)
    level0 = ds.covering_grid(level=0, left_edge=[0,0.0,0.0], dims=ds.domain_dimensions)

    err = np.array([])
    # loop over z slices
    for k in range(ds.domain_dimensions[2]):
        sim = level0['density'][:,:,k].value
        err = np.concatenate([err, (sim - rho0).flatten()])
    rmse.append(np.sqrt(np.mean(err**2)))

outfile = os.path.join(results_dir, 'rmse.csv')
itidx = pd.Index(it,name='iteration')
pd.DataFrame({'t':t,'RMSE':rmse},index=itidx).to_csv(outfile)
print('wrote',outfile)
    
# plot
fig,ax = plt.subplots()
#ax.semilogy(t, rmse, linestyle='none', marker='.')
#ax.set_xlabel('time [s]')
ax.semilogy(it, rmse, linestyle='none', marker='.')
ax.set_xlabel('iteration',fontsize='x-large')
ax.set_ylabel('density RMSE',fontsize='x-large')
outfile = os.path.join(results_dir, 'rmse.png')
fig.savefig(outfile,bbox_inches='tight')
print('saved',outfile)

