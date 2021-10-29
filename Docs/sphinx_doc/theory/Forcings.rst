
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

Pressure Gradient
-----------------

If

::

      abl_driver_type == "PressureGradient"

then

.. math::

  \mathbf{F} = (gpx_{ext}, gpy_{ext}, gpz_{ext})

where :math:`(gpx_{ext}, gpy_{ext}, gpz_{ext})` are user-specified.

Geostrophic Forcing
-------------------

If

::

      abl_driver_type == "GeostrophicWind"

then geostrophic forcing is included in the forcing terms, i.e.

.. math::

  \mathbf{F} = (-C_f \; v_{geo}, C_f \; u_{geo}, 0)

where :math:`C_f = 4 \pi / P_{rot}` is the Coriolis factor with :math:`P_{rot}` the rotational
period (measured in seconds), and :math:`(u_{geo}, v_{geo}, 0)` is the
user-specified geostrophic wind.  Note that if geostrophic forcing is enabled,
Coriolis forcing must also be included.

