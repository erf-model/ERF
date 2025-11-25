# ABL Single Column with Local Timestepping

This test case demonstrates the use of local timestepping to accelerate
convergence to steady state in an idealized neutral atmospheric boundary layer (ABL).

## Description

The simulation uses a single vertical column (1x1 in horizontal) with vertical grid
stretching to resolve the boundary layer. The ABL is driven by a fixed geostrophic
wind which is in balance with Coriolis forces at steady state.

### Configuration

**Domain:**
- Single column: 1x1 horizontal cells, 64 vertical cells
- Height: 2000 m
- Vertical grid stretching ratio: 1.03 (finer near surface)
- Initial vertical spacing: 10 m at surface

**Boundary Conditions:**
- Bottom: MOST (Monin-Obukhov Similarity Theory) surface layer with z₀ = 0.1 m
- Top: Slip wall with fixed θ gradient = 0.003 K/m (matching free atmosphere lapse rate)
- Lateral: Periodic (not relevant for single column)

**Initial Conditions:**
The atmosphere is initialized with an input sounding:
- 0-750 m: Constant potential temperature θ = 300 K (neutral layer)
- 750-850 m: Capping inversion with 8 K temperature increase (100 m thick)
- Above 850 m: Free atmosphere with lapse rate of 3 K/km
- Initial wind: 5 m/s in x-direction (below geostrophic wind to spin up)

**Forcing:**
- Geostrophic wind: 5 m/s in x-direction
- Coriolis force with latitude = 45°

**Turbulence:**
- MRF (Medium Range Forecast) PBL scheme
- Small initial temperature perturbations (0.1 K) to trigger development

## Local Timestepping

With `erf.use_local_timestepping = true`, each vertical level uses its own
timestep based on local CFL constraints. This accelerates convergence because:

- Near-surface cells with high shear advance with smaller timesteps
- Cells in the free atmosphere with lower velocities can use larger timesteps
- The system reaches steady state faster than with a global minimum timestep

## Expected Behavior

At steady state, the flow should reach a balance between:
- Geostrophic wind forcing
- Coriolis force
- Surface friction (via MOST)
- Turbulent mixing (via MRF PBL)

The velocity profile should show:
- Surface layer with logarithmic velocity profile
- Mixed layer with relatively uniform wind
- Transition to geostrophic wind above the boundary layer

## Usage

Run this test case to see local timestepping accelerate convergence to
the equilibrium ABL state. Compare runtime with and without local timestepping
by toggling `erf.use_local_timestepping`.
