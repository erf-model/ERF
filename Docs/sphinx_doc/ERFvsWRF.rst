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
acoustic and gravity wave modes.  Variable time step capability.

**Spatial Discretization**: 2nd- to 6th-order advection options in horizontal and vertical.  In addition, several
different WENO schemes are available for scalar variables other than density and potential temperature.

**Turbulent Mixing**: Sub-grid scale turbulence formulation.  Vertically implicit acoustic step off-centering.

**Diffusion**: In WRF, the diffusion coefficients specified in the input file (:math:`K_h` and :math:`K_v` for
horizontal and vertical diffusion) get divided by the Prandtl number for the potential temperature and the scalars.
For the momentum, they are used as it is. In ERF, the coefficients specified in the inputs (:math:`\alpha_T` and :math:`\alpha_C`)
are used as it is, and no division by Prandtl number is done.

**Initial conditions**: both ERF and WRF have the ability to initialize problems from
3D "real" data (output of real.exe), "ideal" data (output of ideal.exe) and from 1D input soundings.

**Lateral boundary conditions**: Periodic, open, symmetric and specified (in wrfbdy* files).

**Bottom boundary conditions**: Frictional or free-slip

**Earth's Rotation**: Coriolis terms in ERF controlled by run-time input flag

**Mapping to Sphere**: ERF supports the use of map scale factors for isotropic projections (read in from
wrfinput files).

**Nesting**: One-way or two-way.  Multiple levels and integer ratios.



Key Differences
--------------------

**Vertical coordinates**: Unlike WRF, ERF uses a terrain-following height-based vertical coordinate,
with vertical grid stretching permitted.

**Time Integration**: ERF supports using a 3rd-order Runge-Kutta scheme with no substepping as alternative to RK3 with acoustic substepping.

**Initial conditions**: ERF has an additional mode of "custom" initialization in which
the user writes the initialization routine.

ERF does *not* have the capability for global simulation

