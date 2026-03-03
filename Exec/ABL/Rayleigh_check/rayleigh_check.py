#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

ubar = 10.
Rw = 0.02

profiles = np.loadtxt('mean_profiles.dat')
level = profiles[0::4]
t = level[:,0]
u = level[:,2]

fig,ax = plt.subplots()
ax.plot(t,u,'.-',label='simulation')
ax.plot(t,ubar*(1-np.exp(-Rw*t)), 'k--', label=r'analytical solution')
ax.legend()
ax.set_xlabel('t [s]')
ax.set_ylabel('u [m/s]')

plt.show()
