Propagation Methods
===================

The fire front propagation method is selected via::

    erf.fire.propagation_method = farsite   # default
    erf.fire.propagation_method = levelset

FARSITE Ellipse (default)
-------------------------

Lagrangian point-stamping scheme using Richards (1990) elliptical spread
geometry with Anderson (1983) L/W ratio. Points are stamped at the fire front
boundary and advanced according to the local rate of spread, creating a
discrete representation of the fire perimeter.

**Characteristics:**
- Lagrangian point-based tracking of fire front
- Deterministic and stable
- Backwards compatible with all existing input files
- Controlled by ``erf.fire.farsite.*`` parameters

**Parameters:**

- ``erf.fire.farsite.cfl_fire`` (default: 0.5)
  — CFL number for fire subcycle timestep control

- ``erf.fire.farsite.coeff_a``, ``coeff_b``, ``coeff_c`` (default: 0.5, 0.25, 0.1)
  — Richards ellipse coefficients (automatically derived from Anderson L/W if ``use_anderson_lw=1``)

- ``erf.fire.farsite.gaussian_sigma`` (default: −1, auto-derived)
  — Gaussian stamping radius [m]

- ``erf.fire.farsite.phi_threshold`` (default: 0.1)
  — Level-set contour value for front detection


Level-Set Advection (WENO5-Z + SSP-RK3 + Sussman Reinitialization)
-------------------------------------------------------------------

PDE-based propagation that solves the level-set advection equation:

.. math::

    \frac{\partial \phi}{\partial t} = -R(\mathbf{x}) \left( |\nabla \phi| - \varepsilon \Delta \phi \right)

where φ is the signed-distance function (negative inside fire, positive outside),
R(x) is the local rate of spread, ε is an artificial viscosity coefficient,
and Δφ is the Laplacian.

**Numerical scheme:**
- **Spatial:** WENO5-Z (5th-order Weighted Essentially Non-Oscillatory with Z-weighting)
  for reconstruction of flux derivatives
- **Temporal:** Strong Stability Preserving RK3 (3-stage SSP-RK3) for time integration
- **Reinitialization:** Sussman's pseudo-time method to maintain signed-distance property

**Characteristics:**
- Eulerian grid-based implicit level-set tracking
- Topologically clean signed-distance function (φ ≈ ±1 away from front)
- Better mass conservation of fire area compared to Lagrangian stamping
- Suitable for AMR refinement and curvature-dependent models
- Higher computational cost per timestep than FARSITE

**Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``erf.fire.levelset.cfl``
     - 0.4
     - CFL number for fire subcycle timestep (dt_ls = cfl × min(dx,dy) / max_ROS)
   * - ``erf.fire.levelset.eps_visc``
     - 0.4
     - Artificial viscosity coefficient ε in the level-set RHS
   * - ``erf.fire.levelset.reinit_every``
     - 5
     - Reinitialization period: reinitialize φ every N fire subcycles
   * - ``erf.fire.levelset.reinit_iters``
     - 10
     - Number of Sussman reinitialization pseudo-time iterations per reinit
   * - ``erf.fire.levelset.reinit_dtau``
     - −1 (auto)
     - Pseudo-timestep for reinitialization; negative means auto = 0.5·min(dx, dy)


When to use which method
------------------------

**Use FARSITE (default) if:**
- You have existing simulations or benchmarks with FARSITE
- You want maximum stability and predictability
- Computational cost is a primary concern
- The application does not require a strict signed-distance function
- You are using point tracking for visualization or diagnostics

**Use Level-Set if:**
- A topologically clean signed-distance function is required
  (e.g., for curvature-dependent models, signed-distance-based diagnostics)
- Better fire-front area conservation is desired
- You are using AMR mesh refinement and need to resolve the fire at multiple levels
- You prefer Eulerian implicit tracking over Lagrangian point stamping
- Computational cost is not a critical constraint

**Interoperability:**
Both methods share the same underlying ROS computation (Rothermel, MacArthur, etc.)
and are compatible with all other ERF-Fire features (spotting, terrain coupling,
heat flux, etc.). You can switch between methods by changing the
``erf.fire.propagation_method`` parameter without modifying the rest of the setup.


References
----------
- Richards, G. D. (1990). "An elliptical growth model of forest fire fronts and its
  numerical solution." International Journal for Numerical Methods in Engineering, 30(6).
- Anderson, H. E. (1983). "Predicting wind-driven wild fires." USDA For. Serv. Gen.
  Tech. Rep. INT-143.
- Sussman, M. (2003). "A second-order coupled level set and volume-of-fluid method
  for computing growth and collapse of vapor bubbles." Journal of Computational
  Physics, 187(1).
- Osher, S., & Sethian, J. A. (1988). "Fronts propagating with curvature-dependent
  speed: Algorithms based on Hamilton-Jacobi formulations." Journal of Computational
  Physics, 79(1).
