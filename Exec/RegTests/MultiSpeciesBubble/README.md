# MultiSpecies Bubble

Bryan and Fritsch (2002) moist rising-bubble benchmark, run with the super-droplet
(`SuperDroplets`) microphysics to exercise its multi-species capability. Each species is
an independent population of super-droplets condensing from its own water-vapour field;
the second species (`agua`) duplicates the physics of the first (`water`) so that the one-
and two-species runs can be compared directly. The domain is 20 x 0.4 x 10 km on a
200 x 4 x 100 grid (effectively two-dimensional).

## Reference

- Bryan, G.H. and Fritsch, J.M. (2002), *A Benchmark Simulation for Moist Nonhydrostatic
  Numerical Models*, Mon. Wea. Rev., **130**, 2917-2928.
  <https://doi.org/10.1175/1520-0493(2002)130%3C2917:ABSFMN%3E2.0.CO;2>

## Input files

- `inputs_moist_bubble_1species` - moist benchmark bubble, single water species.
- `inputs_moist_bubble_2species` - moist benchmark bubble, two species (`water` + `agua`).
- `inputs_dry_bubble_1species` - dry (isentropic) bubble with one water species.
- `inputs_dry_bubble_2species` - dry (isentropic) bubble with two species.
