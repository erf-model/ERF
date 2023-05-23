#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

zhi = 1000.
Nz = 128
U0 = 5.0
dUdz = 10.0/zhi # (m/s)/m
Ck = 0.1
Ce = 0.93
output_int = 1 # for plot

dz = zhi / Nz

e_exact = Ck/Ce*(dUdz*dz)**2

# verify
e = 0.1
dt = 0.05
dedt = np.nan
earr = []
#for i in range(999999):
for i in range(100):
    if (i%output_int==0):
        earr.append(e)
        print(i,i*dt,dedt,earr[-1])
    # Km == 2 * nu_t
    Km = Ck * dz * np.sqrt(e)
    upwp = -Km * dUdz
    diss = Ce * e**1.5 / dz
    dedt = -upwp*dUdz - diss
    e+= dedt * dt
    if np.abs(dedt) < 1e-8:
        break
print(f"step {i} : nu_t={Km/Km/22} u'w'={upwp} P={-upwp*dUdz} ε={diss} de/dt={dedt} e={e}")
print('expected e =',e_exact)
print('final error:',e - e_exact)

plt.plot(output_int*dt*np.arange(len(earr)),earr)
plt.xlabel('time [s]')
plt.ylabel('SFS e [m$^2$ s$^{-2}$]')
plt.savefig('evolve_e.png',bbox_inches='tight')

#
# expected output:
# e = 0.000656292002688172
#
# dt=1, tol=1e-12
# step 6656 : u'w'=-0.00020014216229914165 P=2.0014216229914164e-06 ε=2.001421632988209e-06 de/dt=-9.996792534261922e-15 e=0.0006562920059562526 (tol=1e-12)
# final error: 3.268080576122878e-12
#
# dt=0.05
# step 99 : nu_t=0.045454545454545456 u'w'=-0.0022610990632662593 P=2.2610990632662593e-05 ε=0.002885902120962711 de/dt=-0.0028632911303300484 e=0.08362116551191935
