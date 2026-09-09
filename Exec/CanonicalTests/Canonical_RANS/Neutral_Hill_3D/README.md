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

The mesh has unit aspect ratio on purpose. With dz different from dx
(20 or 80 m at dx = 40 m) this deck, and the same mesh flattened, produce
a divergence of order 1e139 before the very first projection, the wall
distance Poisson solve then diverges, and the run aborts. The failure is
not the aspect ratio as such (the Askervein deck runs at dx/dz of 0.5, 1
and 2), not the lateral boundaries, the terrain source, grid stretching
or the box layout, and not the RANS code (it appears in the initial
projection, before any turbulence call, and with the closure switched to
Smagorinsky). It is deterministic: under `amrex.init_snan = 1` with the
invalid-operation trap armed, initialisation completes without a trap
and the same 1.788e139 divergence, so it is an arithmetic error tied to
dz relative to dx, not a memory read; a second, separate uninitialised
read in the w boundary fill (`ERFPhysBCFunct_w`) then trips the trap in
the first advance. Reproducer: `inputs_hill3d amr.n_cell="64 64 40"
prob.hmax=1e-6 max_step=0 erf.mg_v=2 erf.v=1`. The
wall distance is `erf.wall_dist_type = terrain_height` (exact to 1e-10 on
a flat mesh, 0.01 % mean error on this hill); the Poisson distance is
exercised by the `_Poisson` CTest variant.

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
