# Askervein Hill

This is an example input file for a neutral atmospheric boundary layer (ABL)
based on the simulation setup described in Wagenbrenner et al. 2019,
Atmosphere.

The modeled domain is a 6 km by 6 km square, centered on the hill.
The grid has Δx = Δy = 20 m in the horizontal and an initial vertical spacing
of Δz = 20 m in the vertical direction with grid stretching ratio of 1.08.
Boundaries are Dirichlet inflow on the west and south, with specified power-law
profile from literature, and Neumann outflow on the north and east; the bottom
boundary is a surface layer based on MOST and the top boundary is a slip wall.

Using the anelastic path with FFT hybrid poisson solver, we perform an unsteady
RANS simulation with a one-equation TKE closure (see Axell & Liungman 2001,
Environ. Fluid Mech.) neglecting Coriolis forces. Damping is provided at the
inlets and top boundary (sponge zone and Rayleigh layer, respectively) to
prevent terrain-induced wave reflections.

The solver executable is compiled in Exec.

## Input decks

The folder holds four decks for the same hill, inflow and stretched mesh. All run the anelastic equations.

| Deck | Model | Levels | Run length |
|---|---|---|---|
| `inputs_anel` | k-equation URANS | 1 | 600 s |
| `inputs_anel_two_level` | k-equation URANS, MidPoint | 2 | 600 s |
| `inputs_les` | Deardorff LES, MidPoint | 1 | 60 s |
| `inputs_les_two_level` | Deardorff LES, MidPoint | 2 | 60 s |

The two-level decks refine a 2 km x 2 km patch around the hill to 10 m (`amr.ref_ratio_vect = 2 2 1`). The vertical direction is stretched, so it is not refined, and the patch must span all of it: the box is given in coarse cell indices (`erf.box1.in_box_lo_indices_crse = 100 100 0`, `in_box_hi_indices_crse = 199 199 17`), because a box given in physical coordinates on a stretched mesh stops short of the top. The decks use `erf.wall_dist_type = terrain_height` for the wall distance and print the divergence after every projection (`erf.mg_v = 1`).

## Refined levels and the projection

The four faces of the refined patch cross the slopes around the hill. The projection on the refined level has zero flux through those faces (see "Neumann faces on terrain-fitted meshes" in the documentation). With the operator ERF used before erf-model/ERF#4279, which let flux through them where the terrain slopes, the refined level never converges and drifts.

## Measured results

The runs used 8 MPI ranks (macOS, Release build with FFT). The table shows, for the refined level, the divergence left after each projection above the constant that a projection on a closed region cannot remove (median and worst over all solves), and the time per solve. The rows marked "before #4279" are the same decks run with the operator ERF used before erf-model/ERF#4279, which is no longer selectable (the same as development at 905623352).

| Deck | Level 1: divergence left (median / worst) | Time per level-1 solve | Theta at the end (both levels) |
|---|---|---|---|
| `inputs_les_two_level` | 2.1e-10 / 2.2e-9 | 0.060 s | 298.7-301.4 K |
| `inputs_les_two_level`, before #4279 | 6.7e-6 / 0.22 | 0.38 s | 260.6-370.8 K |
| `inputs_anel_two_level` | 2.7e-10 / 2.2e-9 | 0.041 s | 299.6-300.4 K |
| `inputs_anel_two_level`, before #4279 | 2.1e-7 / 0.72 | 0.29 s | 299.6-300.4 K |

- With the earlier operator the two-level LES drifts within 60 s: theta leaves the neutral 300 K by up to 70 K. The single-level LES and level 0 converge, so the drift comes from the refined level.
- After 600 s of URANS the two operators differ by up to 0.37 m/s in u and 0.95 m/s in v on the refined level (3 % and 5 % relative), in streaks along the flow downstream of the hill and in the lowest 150 m. Theta hardly differs (6e-3 K on the refined level, 299.88-300.14 K in both), because the URANS case is neutral and its divergence error goes into the wind.
- Level 0 converges either way (worst 8e-9 left).
- `inputs_les` on one level: 259 solves at 0.155 s each, worst 8e-9 left; theta stays at 300 K. The edges of this terrain are flat to round-off, so the two operators agree there to 1e-14 in the wind after the first projection.

### With AMReX's multigrid (`erf.terrain_poisson_solver = mlmg`)

The same decks run with `erf.terrain_poisson_solver = mlmg` added on the command line, in builds with or without FFT. On the same 8 ranks:

| Deck | Level 1: divergence left (median / worst) | Time per level-1 solve |
|---|---|---|
| `inputs_les_two_level` | 3.7e-9 / 1.2e-8 | 0.056 s |
| `inputs_anel_two_level` | 6.6e-10 / 1.3e-8 | 0.054 s |

- The fields match the FFT-GMRES runs to 1.2e-4 K in theta (LES, 60 s) and to 6.6e-5 m/s in the wind and 4e-5 K in theta (URANS, 600 s; theta 299.6-300.4 K in both), 6e-6 relative: the two solvers solve the same equations and differ by their stopping tolerances.
- Level 0: 0.19 s per solve with either solver.
