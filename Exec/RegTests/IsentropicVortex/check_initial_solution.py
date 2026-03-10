#!/usr/bin/env python
import numpy as np
import yt
from isentropic_vortex_helpers import read_input, IsentropicVortex

prob_parm = read_input('inputs_ex',prefix='prob.')
print('Problem parameters from inputs_ex:\n',prob_parm)

# load plotfile
ds = yt.load('plt00000')

# calculate corresponding isentropic vortex solution
vort = IsentropicVortex(**prob_parm)
xx,yy,rho_ref,u_ref,v_ref,p_ref,T_ref = vort.evaluate(
    ds.domain_left_edge,
    ds.domain_right_edge,
    ds.domain_dimensions
)

# load solution at level0 (the only level w/o AMR)
level0 = ds.covering_grid(level=0, left_edge=[0,0.0,0.0], dims=ds.domain_dimensions)

# compare 2D solution on z slices
for k in range(ds.domain_dimensions[2]):
    rmserr = []
    for field,ref in zip(['density','x_velocity','y_velocity','pressure','temp'],
                         [rho_ref,u_ref,v_ref,p_ref,T_ref]):
        sim = level0[field][:,:,k].value
        err = (sim - ref).flatten()
        rmserr.append(np.sqrt(np.mean(err**2)))
    print('k=',k,rmserr)
