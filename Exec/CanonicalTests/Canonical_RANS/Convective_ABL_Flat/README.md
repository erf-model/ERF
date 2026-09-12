# Canonical RANS: convective ABL on flat ground

A 937 m deep neutral layer at 300 K capped by an 8 K inversion (the
Moeng & Sullivan 1994 case B sounding used by `Canonical_LES`), heated
from below by a constant kinematic heat flux of 0.24 K m/s under a 10 m/s
geostrophic wind at f = 1e-4 1/s, z0 = 0.16 m. Closed with the
one-equation k RANS model of Axell & Liungman (2001) under a MOST surface
layer, with `erf.rans_consistent_diffusivities` and the geometric length
capped from the diagnosed boundary-layer height
(`erf.rans_lscale_from_pblh` with `erf.max_geom_lscale = 100` as the
ceiling, since zi is about 1 km).

| item | value |
| --- | --- |
| domain | 2560 x 2560 x 2000 m |
| grid | 8 x 8 x 100, dz = 20 m, first cell centre 10 m |
| time step | 2 s fixed, anelastic with FFT (K/rho reaches 40 m^2/s, the explicit limit at 5 s) |
| closure | `erf.rans_type = kEqn`, AL01 defaults, `dirichlet_k = true` |
| physics run | 4 h (7200 steps), about 2 min on 2 ranks |
| smoke run | 40 steps (`ctest -R RANS_Convective_ABL_Flat`) |

## Running

```bash
mpirun -np 2 erf_exec inputs_convective
python3 check_convective.py --physics plt07200 surf_hist.dat
```

## Checks

Smoke (`--smoke`, the CTest entry): the structural checks of
`../rans_checks.py`, with the length-scale bound taken as the unstable
bound (about 1.31 times the neutral geometric length under the cap), which
is where the phase 3 limiter on the unstable length is exercised.

Physics (`--physics`, 4 h):

| check | target | tolerance |
| --- | --- | --- |
| u* | 0.30 to 0.80 m/s | range |
| column heat gain, sum rho (theta - theta_init) dz, over rho_sfc F t | 1 | 10 % |
| inversion height (strongest dtheta/dz) | 900 to 1250 m | range |
| theta spread between 0.2 zi and 0.7 zi | 0 | 2 K (local-K closure, see below) |
| max dtheta/dz between 0.2 zi and 0.7 zi | <= 0 (superadiabatic or neutral) | 1e-3 K/m |
| mixed-layer warming over the encroachment estimate F t / zi | 1 | 30 % |
| KE at the first cell over u*^2 | >= 1 / Cmu0^2 (buoyancy adds to Eq. 16) | 1 % |
| min KE between 0.1 zi and 0.8 zi | >= 0.05 m^2/s^2 | exact |
| max KE above 1.3 zi | <= 0.05 m^2/s^2 | exact |

The heat-budget check is exact up to the top boundary and the damping
layer; the mixed-layer warming after 4 h is about 3.5 K, which leaves the
inversion at its initial height, so the case tests mixing and the wall
condition rather than entrainment. A local-K closure cannot produce the
countergradient transport of a convective layer: the heat flux F(z) is
carried by K_h dtheta/dz alone, so the profile keeps a superadiabatic
lapse of about -F/K_h, near -2 K/km for K_h of 35 m^2/s, and about 1 K of
spread across the mixed layer where LES shows under 0.3 K. The spread
tolerance bounds that value rather than the LES one.
