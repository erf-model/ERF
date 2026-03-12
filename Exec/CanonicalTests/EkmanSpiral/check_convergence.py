#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yt

v0 = 15.0 # geostrophic wind, from prob.V_0 or erf.geostrophic_wind
Az = 5.0 # kinematic viscosity [m^2/s], from erf.dynamicViscosity / prob.rho_0
earth_period = 86164.0900027328 # erf.rotational_time_period
ymax_plot = None #2000.0

coriolis_factor = 4 * np.pi / earth_period # [1/s]
DE = np.sqrt(2.0 * Az / coriolis_factor) # Ekman depth [m]
print('Ekman depth =',DE,'m')
def UExact(z):
    return v0 * (1.0 - np.exp(-z/DE) * np.cos(-z/DE))
def VExact(z):
    return -v0 * np.exp(-z/DE) * np.sin(-z/DE)

def load_soln(pltfile):
    print('Loading',pltfile)
    ds = yt.load(pltfile)
    soln = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    zp = np.linspace(ds.domain_left_edge[2],
                     ds.domain_right_edge[2],
                     ds.domain_dimensions[2]+1)
    z = (zp[1:] + zp[:-1]).value / 2
    u = soln['x_velocity'][0,0,:].value
    v = soln['y_velocity'][0,0,:].value
    uerr = u - UExact(z)
    verr = v - VExact(z)
    return pd.DataFrame({'u':u,'v':v,'uerr':uerr,'verr':verr},
                        index=pd.Index(z,name='height'))

#
# read all test cases
#
#pltfile = 'plt00000'
pltfile = 'plt00010'
results = {
    25.0: load_soln(f'convergence/dz25.0/plt00010'),
    12.5: load_soln(f'convergence/dz12.5/plt00020'),
    6.25: load_soln(f'convergence/dz6.25/plt00040'),
    3.125: load_soln(f'convergence/dz3.125/plt00080'),
    1.5625: load_soln(f'convergence/dz1.5625/plt00160'),
    0.78125: load_soln(f'convergence/dz0.78125/plt00320'),
}

#
# plot velocity profiles
#
fig,ax = plt.subplots(ncols=2,sharey=True,figsize=(8,6))
for name,df in results.items():
    ax[0].plot(df['u'],df.index,label=f'dz = {name:g}')
    ax[1].plot(df['v'],df.index,label=f'dz = {name:g}')

ax[0].plot(UExact(df.index),df.index,color='0.5',lw=5,alpha=0.3,label='analytical')
ax[1].plot(VExact(df.index),df.index,color='0.5',lw=5,alpha=0.3,label='analytical')

ax[0].set_ylim((0,ymax_plot))
ax[0].set_xlabel('x-velocity [m/s]', fontsize='x-large')
ax[1].set_xlabel('y-velocity [m/s]', fontsize='x-large')
ax[0].set_ylabel('height [m]',fontsize='x-large')
ax[1].legend(loc='upper right',fontsize='large')
for axi in ax:
    axi.tick_params(labelsize='large')
    axi.grid()

fig.savefig('convergence/velocity_profiles.png',bbox_inches='tight')

ax[0].set_ylim((0,50))
fig.savefig('convergence/velocity_profiles_zoom.png',bbox_inches='tight')
#
# plot error profiles
#
fig,ax = plt.subplots(ncols=2,sharey=True,figsize=(8,6))
for name,df in results.items():
    ax[0].semilogx(df['uerr'].abs(),df.index,label=f'dz = {name:g}')
    ax[1].semilogx(df['verr'].abs(),df.index,label=f'dz = {name:g}')

ax[0].set_ylim((0,ymax_plot))
ax[0].set_xlabel('absolute x-velocity error [m/s]', fontsize='x-large')
ax[1].set_xlabel('absolute y-velocity error [m/s]', fontsize='x-large')
ax[0].set_ylabel('height [m]',fontsize='x-large')
ax[1].legend(loc='upper right',fontsize='large')
for axi in ax:
    axi.tick_params(labelsize='large')
    axi.grid()

fig.savefig('convergence/error_profiles.png',bbox_inches='tight')

#
# plot L2 error norm
#
errors = pd.DataFrame({
    name: np.sqrt(np.sum(df[['uerr','verr']]**2))
    for name,df in results.items()
}).transpose()
print(errors)

fig,ax = plt.subplots()
ax.loglog(errors.index, errors['uerr'], '-o', label='u-error')
ax.loglog(errors.index, errors['verr'], '-o', label='v-error')

y1 = errors['uerr'].iloc[0]
y2 = errors['uerr'].iloc[0]/1024.
x = np.array([0.78125, 25.0])
y = np.array([y2,y1])
ax.loglog(x,y, '-', label='second-order')

ax.legend(loc='upper left',fontsize='large')
ax.set_xlabel(r'$\Delta z$ [m]',fontsize='x-large')
ax.set_ylabel(r'$L_2$ error [m/s]',fontsize='x-large')
ax.tick_params(labelsize='large')
ax.grid()

fig.savefig('convergence/error_vs_dz.png',bbox_inches='tight')

