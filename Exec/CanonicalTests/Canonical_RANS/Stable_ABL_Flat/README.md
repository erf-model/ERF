# Canonical RANS: stable ABL on flat ground (GABLS1)

The GABLS1 intercomparison set-up of Beare et al. (2006, Boundary-Layer
Meteorol. 118, 247-272): an 8 m/s geostrophic wind at 73 N over a
z0 = 0.1 m surface cooled at 0.25 K/h from 265 K, starting from a neutral
265 K layer 100 m deep capped by 0.01 K/m. Closed with the one-equation
k RANS model of Axell & Liungman (2001) under a MOST surface layer, with
the horizontal heat and scalar diffusivities following the AL01 scalar
stability function (`erf.rans_consistent_diffusivities`) and the geometric
length capped from the diagnosed boundary-layer height
(`erf.rans_lscale_from_pblh`, `erf.most.pblh_calc = MYNN25`).

| item | value |
| --- | --- |
| domain | 400 x 400 x 400 m |
| grid | 8 x 8 x 100, dz = 4 m, first cell centre 2 m |
| time step | 2 s fixed, anelastic with FFT |
| closure | `erf.rans_type = kEqn`, AL01 defaults, `dirichlet_k = true` |
| physics run | 9 h (16200 steps), about 3 min on 2 ranks |
| smoke run | 40 steps (`ctest -R RANS_Stable_ABL_Flat`) |

## Running

```bash
mpirun -np 2 erf_exec inputs_stable
python3 check_stable.py --physics plt16200 surf_hist.dat
```

## Checks

Smoke (`--smoke`, the CTest entry): the structural checks of
`../rans_checks.py` (finite fields; KE, Kmv, diss non-negative; wall
distance equal to the cell-centre height; Lturb at most the unstable bound
of the neutral geometric length; dissipation consistent with AL01 Eq. 19 in
interior cells within 5 %; wall-cell k retention within 1 %).

Physics (`--physics`, 9 h):

| check | target | tolerance |
| --- | --- | --- |
| u* | 0.20 to 0.35 m/s | range |
| theta at the first cell minus the imposed surface theta | 0 to 1.5 K | range |
| min dtheta/dz below 200 m | >= 0 | 1e-3 K/m |
| max wind speed over Ug (low-level jet) | >= 1.02 | exact |
| height of the wind maximum | 50 to 300 m | range |
| boundary-layer depth from KE (5 % of the wall value) | 80 to 300 m | range |
| KE at the first cell over u*^2 | 1 / Cmu0^2 = 3.23 | 5 % |
| Lturb over the neutral geometric length, 120 to 300 m | <= 1 | 1e-6 |

The GABLS1 LES ensemble gives u* of 0.26 to 0.30 m/s, a jet of about 9.5
m/s near 150 to 200 m and a boundary layer 150 to 200 m deep after 9 h;
the ranges above bracket what a one-equation RANS on a 4 m grid should
reproduce.
