# ABLWarmBubble

Two start-up options of the ABL problem that WRF's idealized fire cases use,
checked against closed forms by `check_warm_bubble.py`:

- a warm bubble, `prob.bubble_dT`, `prob.bubble_x_c/y_c/z_c`,
  `prob.bubble_x_r/y_r/z_r` (WRF `delt_perturbation`, `x/y/zrad_perturbation`,
  `hght_perturbation`): theta' = bubble_dT cos^2(pi r / 2) at each cell's own
  height, made at constant pressure in compressible runs and at the fixed density
  in anelastic runs;
- a `bulk_coeff` surface layer with an imposed kinematic heat flux
  `erf.most.surf_temp_flux` and no `erf.most.Ch` (WRF `isfflx = 0` with
  `tke_heat_flux`).

The mesh is stretched in z (4 m at the ground, ratio 1.15, 12 levels), so a
bubble placed on a uniform (k + 1/2) dz would miss. The run is calm, has no
diffusion inside the domain and lasts 0.5 s (10 steps), so the bubble barely
moves and the surface heat stays in the bottom cell.

| Deck | What it runs |
|------|--------------|
| `inputs_base` | shared settings; on its own it sets the bubble centre and radii without `prob.bubble_dT` |
| `inputs` | compressible run with the bubble and the imposed flux |
| `inputs_anelastic` | the same with `erf.anelastic = 1` (FFT Poisson solve; registered only in an FFT build) |
| `inputs_keys_without_dT` | start-up check: must stop because `prob.bubble_dT` is not set |

`check_warm_bubble.py <plotfile>` checks theta' at every cell above the bottom
one, the density (constant-pressure deficit, or the base state when anelastic),
the pressure perturbation in the bubble (compressible), and the heat budget
sum(rho theta dV) gain = sum(rho_1 F dA) t. `--self-test` proves that it fails a
fixed-density compressible bubble, a uniform-dz bubble and a lost surface flux.

Measured at t = 0.5 s (1 and 2 ranks give the same numbers):

| Check | `inputs` | `inputs_anelastic` | Tolerance |
|-------|----------|--------------------|-----------|
| max \|theta - 300 - theta'\| [K] | 2.6e-5 | 1.5e-5 | 2e-3 |
| max relative density error | 5.3e-6 | 0 | 2e-5 |
| heat gained / (rho F A t) - 1 | 5.0e-5 | 1.5e-7 | 0.02 |
| max \|p'\| in the bubble [Pa] | 0.48 | (not checked) | 5 |

The previous form of the bubble (rho*theta raised at fixed density in a
compressible run) gives 14 Pa, a 1.7e-4 density error and a heat budget off by
3.2; a bubble on uniform-dz heights misses theta' by 0.48 K; a surface layer that
ignored the imposed flux (Ch form with Ch = 0) adds no heat (budget off by 1).
