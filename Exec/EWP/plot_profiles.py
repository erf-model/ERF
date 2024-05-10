#!/usr/bin/env python
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import yt

pltfiles = [pltfile for pltfile in glob.glob('plt*')
            if os.path.isdir(pltfile) and not ('old' in pltfile)]
nsteps = [int(pltfile[3:]) for pltfile in pltfiles]
latestoutput = pltfiles[np.argmax(nsteps)]
print(latestoutput)

ds = yt.load(latestoutput)

soln = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
zp = np.linspace(ds.domain_left_edge[2],
                 ds.domain_right_edge[2],
                 ds.domain_dimensions[2]+1)
z = (zp[1:] + zp[:-1]).value / 2

utot = soln['x_velocity'].value
vtot = soln['y_velocity'].value
wtot = soln['z_velocity'].value

U = np.mean(utot,axis=(0,1))
V = np.mean(vtot,axis=(0,1))
W = np.mean(wtot,axis=(0,1))
u = utot - U[np.newaxis,np.newaxis,:]
v = vtot - V[np.newaxis,np.newaxis,:]
w = wtot - W[np.newaxis,np.newaxis,:]

uu = np.var(u, axis=(0,1))
vv = np.var(v, axis=(0,1))
ww = np.var(w, axis=(0,1))

uw = np.mean(u*w, axis=(0,1))
vw = np.mean(v*w, axis=(0,1))
uv = np.mean(u*v, axis=(0,1))

fig,ax = plt.subplots(nrows=3,ncols=3,sharey=True,figsize=(8.5,11))
ax[0,0].plot(U,z)
ax[0,1].plot(V,z)
ax[0,2].plot(W,z)
ax[0,0].set_xlabel(r'$\langle U \rangle$ [m/s]')
ax[0,1].set_xlabel(r'$\langle V \rangle$ [m/s]')
ax[0,2].set_xlabel(r'$\langle W \rangle$ [m/s]')
ax[1,0].plot(uu,z)
ax[1,1].plot(vv,z)
ax[1,2].plot(ww,z)
ax[1,0].set_xlabel(r"$\langle u'u' \rangle$ [m/s]")
ax[1,1].set_xlabel(r"$\langle v'v' \rangle$ [m/s]")
ax[1,2].set_xlabel(r"$\langle w'w' \rangle$ [m/s]")
ax[2,0].plot(uw,z)
ax[2,1].plot(vw,z)
ax[2,2].plot(uv,z)
ax[2,0].set_xlabel(r"$\langle u'w' \rangle$ [m/s]")
ax[2,1].set_xlabel(r"$\langle v'w' \rangle$ [m/s]")
ax[2,2].set_xlabel(r"$\langle u'v' \rangle$ [m/s]")
ax[0,0].set_ylabel('$z$ [m]')
ax[1,0].set_ylabel('$z$ [m]')
ax[2,0].set_ylabel('$z$ [m]')
ax[0,0].set_ylim((ds.domain_left_edge[2], ds.domain_right_edge[2]))
for axi in ax.ravel():
    axi.grid()

fig.savefig('profiles.png',bbox_inches='tight')

plt.show()

