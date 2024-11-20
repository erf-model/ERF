
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Forcings:

Physical Forcings
=================

Physical forcings available in ERF comprise the standard source terms for atmospheric modeling.
These include Coriolis and geostrophic forcing; Rayleigh damping and sponge layer(s); subsidence;
simplified radiative thermal sources; and solution nudging towards a prescribed input sounding.

ERF also supports models for wind farm parametrization in which the effects of wind turbines are represented
by imposing a momentum sink on the mean flow and/or turbulent kinetic energy (TKE).
Currently the Fitch model, Explicit Wake Parametrization (EWP) model, Simplified Actuator Disk model (SAD),
and Generalized Actuator Disk model (GAD) are supported. See :ref:`sec:WindFarmModels` for more information.

Below is more detail on how to set the forcing terms.

Buoyancy
--------

If

::

      use_gravity == true

then buoyancy is included in the momentum equations.  See :ref:`sec:Buoyancy` for more detail
about the possible formulations of the buoyancy term.

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

Values for ``erf.rotational_time_period``, ``erf.latitude``, and ``erf.coriolis_3d``; the first two are used
to compute the Coriolis frequency and the last of these determines whether to include the z-component in the Coriolis forcing.

There is no dependence on the radial distance from the center of the earth, thus the curvature of the earth is neglected.

Rayleigh Damping
----------------

Rayleigh damping can be imposed on any or all of :math:`u, v, w, T` and is controlled by
setting

::

      rayleigh_damp_U = true
      rayleigh_damp_V = true
      rayleigh_damp_W = true
      rayleigh_damp_T = true

in the inputs file.  When one or more of those is true,
explicit Rayleigh damping is included in the energy and/or momentum equations
as described in Section 4.4.3 of the WRF Model Version 4 documentation (p40), i.e. :

.. math::

  \mathbf{F} = - \tau(z) \rho \; (u - \overline{u}, v - \overline{v}, 0)

and

.. math::

  F_{\rho \theta} = - \tau(z) \rho (\theta - \overline{\theta})

where :math:`(\overline{u}, \overline{v}, 0)` is the reference state velocity, typically
defined as the initial horizontally homogeneous fields in idealized simulations,
and :math:`\overline{\theta}` is the reference state potential temperature.
As in the WRF model, the reference state vertical velocity is assumed to be zero.

Sponge regions
----------------------

ERF provides the capability to apply sponge source terms near domain boundaries to prevent spurious reflections that otherwise occur
at the domain boundaries if standard extrapolation boundary condition is used. The sponge zone is implemented as a source term
in the governing equations, which are active in a volumteric region at the boundaries that is specified by the user in the inputs file.
Currently the target condition to which the sponge zones should be forced towards is to be specified by the user in the inputs file.

.. math::

   \frac{dQ}{dt} = \mathrm{RHS} - A\xi^n(Q-Q_\mathrm{target})

where RHS are the other right-hand side terms. The parameters to be set by the user are -- `A` is the sponge amplitude, `n` is the sponge strength and the :math:`Q_\mathrm{target}` -- the target solution in the sponge. :math:`\xi` is a linear coordinate that is 0 at the beginning of the sponge and 1 at the end. An example of the sponge inputs can be found in ``Exec/RegTests/Terrain2d_Cylinder`` and is given below. This list of inputs specifies sponge zones in the inlet and outlet of the domain in the x-direction and the outlet of the domain in the z-direction. The `start` and `end` parameters specify the starting and ending of the sponge zones. At the inlet, the sponge starts at :math:`x=0` and at the outlet the sponge ends at :math:`x=L` -- the end of the domain. The sponge amplitude `A` has to be adjust
ed in a problem-specific manner. The density and the :math:`x, y, z` velocities to be used in the sponge zones have to be specified in the inputs list.

::

          erf.sponge_strength = 10000.0
          erf.use_xlo_sponge_damping = true
          erf.xlo_sponge_end = 4.0
          erf.use_xhi_sponge_damping = true
          erf.xhi_sponge_start = 26.0
          erf.use_zhi_sponge_damping = true
          erf.zhi_sponge_start = 8.0

          erf.sponge_density = 1.2
          erf.sponge_x_velocity = 10.0
          erf.sponge_y_velocity = 0.0
          erf.sponge_z_velocity = 0.0

Another way of specifying sponge zones is by providing the sponge zone data as a text file input. This is currently implemented only for forcing :math:`x` and :math:`y` velocities in the sponge zones.
The sponge data is input as a text file with 3 columns containing :math:`z, u, v` values. An example can be found in ``Exec/SpongeTest`` and a sample inputs list for using this feature is given below. This list specifies a sponge zone in the inlet in the x-direction. The :math:`u` and :math:`v` velocity forcing in the sponge zones will be read in from the text file -- `input_sponge_file.txt`.

::

          erf.sponge_type = "input_sponge"
          erf.input_sponge_file = "input_sponge_file.txt"
          erf.sponge_strength = 1000.0
          erf.use_xlo_sponge_damping = true
          erf.xlo_sponge_end = 4.0

Problem-Specific Forcing
========================

The following two options can be used to specify external forcing terms.

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

