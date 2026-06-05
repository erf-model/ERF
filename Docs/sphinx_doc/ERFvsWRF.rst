 .. role:: cpp(code)
    :language: c++

.. _ERFvsWRF:

ERF vs WRF
===============

The following comparison is based on the WRF Version 4 Technical Report, titled
"A Description of the Advancd Research WRF Model Version 4"

Similarities
--------------------

**Equations**: both ERF and WRF solve the fully-compressible, Eulerian nonhydrostatic equations, and conserve
dry air mass and scalar mass.  ERF does not have a hydrostatic option.

**Prognostic Variables**: velocity components (u,v,w); perturbation moist potential temperature.  Optionally,
turbulent kinetic energy and any number of scalars such as water vapor mixing ratio, rain/snow mixing ratio,
cloud water / ice mixing ratio.

**Horizontal grid**: both ERF and WRF use Arakawa C-grid staggering.

**Time Integration**: Time-split integration using 3rd-order Runge-Kutta scheme with smaller time step for
acoustic and gravity wave modes.  Variable time step capability. Vertically implicit acoustic step off-centering.

**Spatial Discretization**: 2nd- to 6th-order advection options in horizontal and vertical.  In addition, several
different WENO schemes are available for scalar variables other than density and potential temperature.

**Turbulent Mixing**: ERF and WRF have the same sub-grid scale turbulence closures with the Smagorinsky or
1.5-order TKE (Deardorff) model, in isotropic or anisotropic forms, for large-eddy simulation (LES);
planetary boundary layer (PBL) schemes (MYNN, YSU) are available. ERF also has support for RANS turbulence modeling.

**Diffusion**: In WRF and ERF, constant diffusion coefficients may be specified (:math:`K_h` and :math:`K_v` for
horizontal and vertical diffusion). Constant dynamic viscosity may also be specified in ERF.
Variable diffusivity is provided in 3-D through LES modeling and in 1-D through PBL modeling. For mesoscale applications,
3-D diffusion is provided by combining a PBL scheme with the Smagorinsky model.
Prandtl and Schmidt numbers are used to derive diffusivities of heat or other scalars from the diffusivity of momentum.

**Initial conditions**: both ERF and WRF have the ability to initialize problems from
3-D "real" data (output of real.exe), "ideal" data (output of ideal.exe), and from 1-D input soundings.

**Lateral boundary conditions**: Periodic, open, symmetric and specified (in wrfbdy* files).

**Bottom boundary conditions**: Frictional (Monin-Obukhov Similarity Theory) or free-slip

**Earth's Rotation**: Coriolis terms in ERF controlled by run-time input flag (2-D or 3-D, constant or spatially varying
for real-data cases)

**Mapping to Sphere**: ERF supports the use of map scale factors for isotropic projections (read in from
wrfinput files).

**Nesting**: One-way or two-way.  Multiple levels and integer ratios.

**Wind Energy Modeling**: Wind farm parameterizations and a generalize actuator disk are available.


Key Differences
--------------------
ERF provides **performance portability** on different computing architectures **including GPUs from all major vendors** (NVIDIA, AMD, and Intel).

**Vertical Coordinates**: Unlike WRF, ERF uses a height-based vertical coordinate, with vertical grid stretching permitted.

**Governing Equations**: ERF supports both fully compressible and anelastic equation sets.

**Time Integration**: ERF supports using a 3rd-order Runge-Kutta scheme with explicit acoustic substepping or no substepping
(in addition to the implicit acoustic substepping in WRF).

**Representation of Surface Features**: Terrain and urban geometries may be simulated with immersed forcing or embedded (immersed) boundary techniques,
in addition to the terrain-fitted coordinates approach.

**Interface with AMR-Wind**: ERF may be tightly coupled with AMR-Wind, an incompressible ABL solver with integrated turbine aeroservoelastic dynamics modeling and two-phase flow capabilities.

**Particles**: ERF can be compiled with support for particles, for flow visualization or Lagrangian physics modeling.

**User-Defined Functions**: ERF provides templates to customize initialization and/or impose spatiotemporally varying source terms.

ERF does *not* have the capability for global simulation

