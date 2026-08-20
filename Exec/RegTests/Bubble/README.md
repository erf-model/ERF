# Bubble

A collection of two- and three-dimensional buoyant-bubble and related benchmarks built on
the single `Bubble` problem setup. The cases share the same problem source but differ in
initialization, microphysics, refinement, and boundary conditions. The two-dimensional
cases use a thin (4-cell) y-direction. The dry and moist bubble benchmarks follow Bryan
and Fritsch (2002).

## Reference

- Bryan, G.H. and Fritsch, J.M. (2002), *A Benchmark Simulation for Moist Nonhydrostatic
  Numerical Models*, Mon. Wea. Rev., **130**, 2917-2928.
  <https://doi.org/10.1175/1520-0493(2002)130%3C2917:ABSFMN%3E2.0.CO;2>

## Dry bubble (Bryan-Fritsch 2002)

- `inputs_BF02_dry_bubble` - 2D dry isentropic rising bubble (20 x 10 km, 200 x 100), no
  microphysics; the reference dry benchmark.
- `inputs_BF02_dry_bubble_AMR1`, `inputs_BF02_dry_bubble_AMR2` - same 2D bubble with one
  and two levels of mesh refinement and super-droplet particles, exercising SDM particle
  advection across refinement boundaries.
- `inputs_3DBF02_dry_bubble` - 3D dry bubble (20 x 20 x 9.6 km, 100 x 100 x 48) with
  super-droplet particles.
- `inputs_3DBF02_dry_bubble_AMR1`, `inputs_3DBF02_dry_bubble_AMR2` - 3D dry bubble with one
  and two refinement levels.

## Moist bubble (Bryan-Fritsch 2002)

- `inputs_BF02_moist_bubble` - 2D moist benchmark bubble with the bulk Kessler scheme.
- `inputs_BF02_moist_bubble_SDM` - same moist bubble with super-droplet microphysics
  (NaCl aerosol).
- `inputs_BF02_moist_bubble_SDM_multi_injections_unimodal_NaCl` - moist SDM bubble with
  three particle-injection sources (two moving box regions with opposing velocities and one
  time-limited bubble injection); a demonstration of injection and moving source regions.

## Other cases

- `inputs_grav2d_x` - two-dimensional gravity current driven by a buoyancy perturbation
  (51.2 x 6.4 km, 512 x 64).
- `inputs_squall2d_x` - two-dimensional squall line (50 x 20 km, 201 x 80) with the SAM
  scheme, initialized from `input_sounding_squall2d`.
- `inputs_test_outflow` - outflow boundary-condition test on a moist Kessler bubble.

## Sounding file

- `input_sounding_squall2d` - sounding for `inputs_squall2d_x`.
