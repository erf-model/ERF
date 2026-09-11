# Canonical RANS: neutral flow over an axisymmetric hill

The three-dimensional counterpart of `Neutral_Hill_2D`: a periodic radial
Witch of Agnesi hill, hmax = 100 m and L = 500 m (half-width at half
height), on a terrain-fitted mesh under a MOST surface layer, driven by a
10 m/s geostrophic wind at f = 1e-4 1/s, z0 = 0.1 m, closed with the
one-equation k RANS model of Axell & Liungman (2001). The wall distance
is checked against the exact distance to the surface of revolution and
the crest speed-up against the axisymmetric estimate of about 1.6 h/L
(the 2D Jackson-Hunt value 2 h/L reduced by the flow around the sides).

| item | value |
| --- | --- |
| domain | 5120 x 5120 x 800 m |
| grid | 64 x 64 x 20, dx = dy = dz = 40 m, basic terrain following |
| time step | 1.5 s fixed (Courant 0.375), anelastic |
| closure | `erf.rans_type = kEqn`, AL01 defaults, `dirichlet_k = true` |
| physics run | 4 h (9600 steps), about 10 min on 2 ranks |
| smoke run | 40 steps (`ctest -R RANS_Neutral_Hill_3D`) |

The mesh has unit aspect ratio. The 1.788e139 divergence before the
first projection that this deck showed with dz different from dx
(`amr.n_cell = 64 64 40` at `amr.max_grid_size = 32`) was not the aspect
ratio: that mesh splits the BoxArray in z, and the initial projection
read the unfilled momenta ghost faces at the internal box faces while
the planar surface-layer arrays were duplicated across the stacked boxes
(fixed in erf-model/ERF#3970). What does depend on dz relative to dx is
the Poisson wall-distance solve, whose multigrid diverges at dx = 2 dz
with or without the split, so the `_Poisson` CTest variant needs dz = dx;
the `terrain_height` distance the deck uses has no such limit and is
gathered from the surface boxes, so a z-split layout is fine with it.
The wall distance is `erf.wall_dist_type = terrain_height` (exact to
1e-10 on a flat mesh, 0.01 % mean error on this hill); the Poisson
distance is exercised by the `_Poisson` CTest variant.

## Running

```bash
mpirun -np 2 erf_exec inputs_hill3d
python3 check_hill3d.py --physics plt09600
```

## Checks

Wall distance (both modes, on every second column):

| check | target | tolerance |
| --- | --- | --- |
| max relative error vs the exact distance to the hill | 0 | 15 % (crest cells) |
| mean relative error | 0 | 3 % |
| max absolute error where the distance is under 100 m | 0 | 0.2 dz |

Structural (both modes, whole field): finite fields, KE and Kmv
non-negative, Lturb at most the unstable bound of the geometric length
built on the wall distance, wall-cell k retention within 1 % over the
whole surface.

Physics (`--physics`, 4 h):

| check | target | tolerance |
| --- | --- | --- |
| crest speed-up U_top/U_up - 1 in the three lowest cells | 1.6 h/L = 0.32 | 0.16 to 0.64 |
| speed-up positive in the lowest eight cells | > 0 | exact |
| upstream u* from the wall k | 0.25 to 0.55 m/s | range |
| upstream wind in the three lowest cells vs the log law | equal | 15 % |
