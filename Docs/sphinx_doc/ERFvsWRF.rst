 .. role:: cpp(code)
    :language: c++

.. _ERFvsWRF:

ERF vs WRF
===============

The following comparison is based on the WRF Version 4 Technical Report, titled
"A Description of the Advancd Research WRF Model Version 4"

Similarities
--------------------

Equations: both ERF and WRF solve the fully-compressible, Eulerian nonhydrostatic equations, and conserve
dry air mass and scalar mass.  ERF does not have a hydrostatic option.

Prognostic Variables: velocity components (u,v,w); perturbation moist potential temperature.  Optionally,
turbulent kinetic energy and any number of scalars such as water vapor mixing ratio, rain/snow mixing ratio,
cloud water / ice mixing ratio.

Horiontal grid: both ERF and WRF use Arakawa C-grid staggering.

Time Integration: Time-split integration using 3rd-order Runge-Kutta scheme with smaller time step for
acoustic and gravity wave modes.  Variable time step capability.

Spatial Discretization: 2nd- to 6th-order advection options in horizontal and vertical

Turbulent Mixing: Sub-grid scale turbulence formulation.  Vertically implicit acoustic step off-centering.

Initial conditions: both ERF and WRF have the ability to initialize problems from 3D "real" data (output of real.exe),
"ideal" data (output of ideal.exe) and from 1D input soundings.  In addition, ERF simulations can be initialized
in "custom" mode wherein the user writes the initialization routine.

Lateral boundary conditions: Periodic, open, symmetric and specified (in wrfbdy* files).

Bottom boundary conditions: Frictional or free-slip

Earth's Rotation: Coriolis terms in ERF controlled by run-time input flag

Nesting: One-way or two-way.  Multiple levels and integer ratios.


Key Differences
--------------------

Vertical coordinates: Unlike WRF, ERF uses a terrain-following height-based vertical coordinate, with vertical
grid stretching permitted.
