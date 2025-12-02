
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

then buoyancy is included in the momentum equations.  See :ref:`Buoyancy` for more detail
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

When initializing from a ``wrfinput`` or ``met_em`` file, the latitude at the grid cell centers will be known. For this case, a user may specify

::

      variable_coriolis == true

to use the grid latitude, :math:`\phi(y)`, when computing the sine and cosine coefficients above.

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

in the inputs file.  When one or more of those is true and
the Rayleigh damping type is set to SlowExplicit or FastExplict,
then explicit Rayleigh damping is included in the energy and/or momentum equations
in the form described in Section 4.4.3 of the WRF Model Version 4 documentation (p40), i.e. :

.. math::

  \mathbf{F} = - \tau(z) \rho \; (u - \overline{u}, v - \overline{v}, w - 0)

and

.. math::

  F_{\rho \theta} = - \tau(z) \rho (\theta - \overline{\theta})

where :math:`(\overline{u}, \overline{v}, 0)` is the reference state velocity, typically
defined as the initial horizontally homogeneous fields in idealized simulations,
and :math:`\overline{\theta}` is the reference state potential temperature.
As in the WRF model, the reference state vertical velocity is assumed to be zero.

If the Rayleigh damping type is set to SlowExplicit then all the damping terms are computed once per
RK stage; if the type is FastExplicit then the damping terms are computed once per acoustic substep.
Either way, they are added explicitly.

If the Rayleigh damping type is set to FastImplicit then the damping term for w only is included implicitly
within the acoustic substepping algorithm; any additional Rayleigh damping (e.g. for u, v, or T) occurs
as if the type is FastExplicit.

The algorithm for FastExplicit is as described in (3.44) of the `MPAS report`_
which is equivalent to that written in (9) of
`Klemp, Dudhia & Hassiotis, An Upper Gravity-Wave Absorbing Layer for NWP Applications (2008)`_.

.. _`MPAS report`: https://www2.mmm.ucar.edu/projects/mpas/mpas_website_linked_files/MPAS-A_tech_note.pdf

.. _`Klemp, Dudhia & Hassiotis, An Upper Gravity-Wave Absorbing Layer for NWP Applications, 2008`: https://journals.ametsoc.org/view/journals/mwre/136/10/2008mwr2596.1.xml

Sponge regions
----------------------

ERF provides the capability to apply sponge source terms near domain boundaries to prevent spurious reflections that otherwise occur
at the domain boundaries if standard extrapolation boundary condition is used. The sponge zone is implemented as a source term
in the governing equations, which are active in a volumetric region at the boundaries that is specified by the user in the inputs file.
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


Immersed forcing to represent terrain
----------------------

An additional option for representing terrain in ERF is to use an immersed forcing method where large body forces are applied to the momentum equations as sinks to force the velocity to near zero or to a desired value.
This method follows the methods of `Chan and Leach (2007) <https://doi.org/10.1175/2006JAMC1321.1>`_ and `Muñoz-Esparza et al. (2020) <https://doi.org/10.1029/2020MS002141>`_, but is expanded to allow the user to utilize a wall-model (based on Monin Obukhov similarity theory, :ref:`sec:surface_layer`).
During initialization, we determine a mask (:math:`\beta_r`) over the entire domain by calculating each cell's volume fraction, which indicates how much of a cell is filled by terrain.
Fully immersed cells have a value of 1, free cells have a value of 0, and partially immersed cells have a value between 0 and 1.
The goal is to force interior cells to near-zero velocities using the following formulation:

.. math::

    F_{\rho u_i} = -C_{d,m} \beta_r \sqrt[3]{\Delta x \Delta y \Delta z} \rho u_i U

where :math:`C_{d,m}` is a drag coefficient and :math:`U` is the wind speed magnitude.
The drag coefficient can be specified by the user using ``erf.if_Cd_momentum``, which defaults to a value of 10.
A larger drag coefficient results in smaller velocities for immersed cells but may require a smaller timestep due to the stiffness of the force.

For partially immersed cells, the user has the option to specify whether to use MOST: ``erf.if_use_most``.
If the user does not specify MOST, then the equation above will be applied and there is an implicit no-slip boundary condition.
If the user does specify MOST, then the following formulation is applied to partially immersed cells:

.. math::

    F_{\rho u_i} = -C_{d,m} (1 - \beta_r) \sqrt[3]{\Delta x \Delta y \Delta z} \rho |U_s| (u_i - u_{i,target})

where :math:`u_{i,target}` is a value determined through MOST and :math:`|U_s|` is a unit velocity scale.
This formulation essentially forces the velocity at the wall to a value determined by using MOST, but the strength forcing is inversely related to how immersed the cell is.
For cells that are more immersed, there is weaker forcing to the target velocity while for cells that are less immersed, there is stronger forcing to the MOST value.

Temperature forcing is also available to represent the temperature of the 'surface'.
The user can specify either a surface temperature and heating rate, a surface flux, or an Obukhov length.
The temperature forcing is then formulated as follows:

.. math::

    F_{\rho\theta} = -C_{d,s} \beta_r \sqrt[3]{\Delta x \Delta y \Delta z} |U_s| (\rho \theta_{target} - \rho\theta)

The target temperature :math:`\theta_{target}`` is straightforward when using a surface temperature and heating rate; when specifying a surface flux or Obukhov length, the target temperature is determined using MOST.
The following inputs are available when representing terrain using immersed forcing:

::

        erf.if_Cd_scalar               = FLOAT
        erf.if_Cd_momentum             = FLOAT
        erf.if_z0                      = FLOAT
        erf.if_surf_temp_flux          = FLOAT
        erf.if_init_surf_temp          = FLOAT
        erf.if_surf_heating_rate       = FLOAT
        erf.if_Olen                    = FLOAT
        erf.if_use_most                = BOOL
        erf.immersed_forcing_substep   = BOOL

An example of using immersed forcing for a Witch of Agnesi hill is available in ``Exec/ABL/immersed_forcing``.

.. note:: When using fully compressible simulations, it is recommended to apply immersed forcing on the substep for numerical stability.

Immersed forcing to represent buildings
----------------------

The immersed forcing capability can also be used to represent buildings.
Currently, the implementation is similar to the formulation for fully immersed cells for terrain, but is proportional to the volume fraction :math:`V_f` thus implicitly applying a no-slip boundary condition following `Muñoz-Esparza et al. (2020) <https://doi.org/10.1029/2020MS002141>`_.
A more advanced wall-model for buildings will be added in the near future (see `Wise et al. (2025) <https://ams.confex.com/ams/25BLT/meetingapp.cgi/Paper/460715>`_ for additional details and a demonstration).
The momentum forcing is defined as follows:

.. math::

    F_{\rho u_i} = -C_{d,m} \beta_r V_f \sqrt[3]{\Delta x \Delta y \Delta z} \rho u_i U

Temperature forcing for building walls and roofs can similar be specified following the same formulation as immersed forcing for terrain; however, only the option to specify a surface temperature and heating rate is currently available.

Inputs that can be used with immersed forcing for buildings are as follows:

::

        erf.buildings_type             = STRING #ImmersedForcing or None
        erf.buildings_file_name        = STRING
        erf.if_Cd_scalar               = FLOAT
        erf.if_Cd_momentum             = FLOAT
        erf.if_init_surf_temp          = FLOAT
        erf.if_surf_heating_rate       = FLOAT
        erf.immersed_forcing_substep   = BOOL

Immersed forcing for buildings can be used in conjunction with the ``StaticFittedMesh`` terrain option.
However, currently, the user must specify the z-coordinates using ``erf.terrain_z_levels``.
In the future, this requirement will be removed.
Note that the volume fraction is calculated prior to the grid transformation; therefore, building heights when located in steep terrain should be considered approximate.

An example of immersed forcing for a building located on top of a Witch of Agnesi hill is available in ``Exec/ABL/immersed_forcing``.


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

