# Isentropic Vortex with Local Timestepping

This test case demonstrates the use of local timestepping to accelerate
convergence to steady state.

## Description

The test uses the isentropic vortex problem, which creates a smooth
vortical flow field. With local timestepping enabled, each cell in the
computational domain uses its own timestep based on local CFL constraints,
allowing cells in regions of lower flow speed to advance faster in time.

## Local Timestepping

The key difference from the standard `IsentropicVortexStationary` test is
the addition of:

```
erf.use_local_timestepping = true
```

This enables per-cell timestep computation based on:
- Local velocity magnitude
- Local sound speed  
- Cell dimensions

## Usage

This feature should only be used for problems where:
1. The steady-state solution is of interest (not transient behavior)
2. Time-accurate resolution is not required
3. The flow field has regions with varying characteristic speeds

## Expected Behavior

- Each cell advances with its own timestep
- Regions with lower speeds converge faster
- The global dt reported is still the minimum across all cells
- Conservation properties are maintained locally
