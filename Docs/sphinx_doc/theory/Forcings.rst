
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Forcings:

Forcings
========

ERF includes the following forcing terms as options:

Buoyancy
--------

If

::

      use_gravity == true

then buoyancy is included in the momentum equations in the form

.. math::

  (0, 0, -\rho^\prime g)

Coriolis Forcing
----------------

If

::

      use_coriolis == true

then Coriolis forcing is included in the momentum equations, i.e. :

.. math::

  \mathbf{F} = (C_f \; (\rho v \sin{\phi} - \rho w \cos{\phi}), -C_f \; \rho u \sin{\phi}, C_f \; \rho u \cos{\phi})

where :math:`C_f = 4 \pi / P_{rot}` is the Coriolis factor with :math:`P_{rot}` the rotational
period (measured in seconds), and :math:`\phi` the latitude.

There is no dependence on the radial distance from the center of the earth, thus the curvature of the earth is neglected. 

Pressure Gradient
-----------------

If

::

      abl_driver_type == "PressureGradient"

then

.. math::

  \mathbf{F} = (\nabla p_{x,ext}, \nabla p_{y,ext}, \nabla p_{z,ext})

where :math:`(\nabla p_{x,ext}, \nabla p_{y,ext}, \nabla p_{z,ext})` are user-specified through ``erf.abl_pressure_grad``.

Geostrophic Forcing
-------------------

If

::

      abl_driver_type == "GeostrophicWind"

then geostrophic forcing is included in the forcing terms, i.e.

.. math::

  \mathbf{F} = (-C_f \; v_{geo}, C_f \; u_{geo}, 0)

where :math:`C_f = 4 \pi / P_{rot}` is the Coriolis factor with :math:`P_{rot}` the rotational
period (measured in seconds), and the geostrophic wind :math:`(u_{geo}, v_{geo}, 0)` is
user-specified through ``erf.abl_geo_wind``.  Note that if geostrophic forcing is enabled,
Coriolis forcing must also be included.

