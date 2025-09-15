# DNS of Turbulent Flow through a Plane Channel

This is an example input file for DNS of channel flow at a Reynolds number $Re_\tau$ of 395 (based
on friction velocity $u_\tau$ and the channel half height $\delta$) -- see Moser, Kim, and Mansour
1999, _Phys. Fluids_. Reference data were download from [this archive](https://turbulence.oden.utexas.edu/MKM_1999.html).

The driving pressure gradient is determined by enforcing a fixed mass flux, based on a bulk
velocity $u/u_\tau = 17.54$. This was determined by integrating the reference streamwise velocity
profile. The corresponding Reynolds number based on channel height and bulk velocity is 13,900
therefore the flow is expected to be fully turbulent.

The solver executable is compiled in `Exec/DryRegTests/Couette_Poiseuille` and provides options for
initializing to a parabolic velocity profile with divergence-free initial perturbations. Random
perturbations were ineffective (flow remains laminar after simulating 20 nondimensional time units).

## Simulation Setup Details

- The solver is run in incompressible mode without buoyancy (gravity is turned off).

- Solver variables are assumed to be nondimensional (indicated by an overbar), normalized by
  $\delta$ and $u_\tau$ as follows:

  $\bar{t} = t u_\tau/\delta$

  $\bar{x} = x / \delta$

  $\bar{u} \equiv u^+ = u / u_\tau$

  $\bar{p} = p / (\rho_\infty u_\tau^2)$

  $\bar{\nu} = 1 / Re_\tau$

- The nondimensional coordinates are related to wall units:
  $x^+ \equiv (u_\tau/\nu) x = Re_\tau \bar{x}$

- We use z as the wall-normal direction instead of y. Grid stretching is defined by a hyperbolic tangent
  function (Abe, Kawamura, and Matsuo 2001, J. Fluids Eng.)

  ```python
  N = 320
  aval = 0.989
  xi = -1 + 2*np.arange(N+1) / N
  zlevels = 1/aval * np.tanh(xi * np.arctanh(aval)) + 1
  print('erf.terrain_z_levels =',' '.join([str(z) for z in zlevels]))
  ```
