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

The solver executable is compiled in Exec/ABL.
