
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _Equations:



Compressible Navier-Stokes Equations
====================================

ERF advances the following set of equations in DNS mode:

.. math::

  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}),

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + pI) +\rho \mathbf{g} + \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho \theta)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \theta) + \nabla \cdot (\alpha_{T}\ \nabla (\rho \theta)),

  \frac{\partial (\rho C)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} C) + \nabla \cdot (\alpha_{C}\ \nabla (\rho C)),

where :math:`\tau` is the stress tensor and :math:`\mathbf{F}` are the forcing terms described in :ref:`Forcings`.

When run in DNS mode, :math:`\tau = 2 \mu S` where :math:`\mu` is the user-specified dynamic viscosity and
:math:`S` is the strain-rate tensor.  When using the Smagorinsky model, the turbulent viscosity
:math:`K = 2 (C_s \Delta)^2 (\sqrt(2 S S) \rho` is used in place of :math:`2 \mu,` where
:math:`C_s` is the Smagorinsky constant and :math:`\Delta` is the mesh spacing.

Note that :math:`\alpha_{T}` is in general variable for a general compressible flow. However, for low Mach number atmospheric flows it are assumed to be constant.

These equations can be re-written in perturbational form by replacing the z-momentum equation with

.. math::

  \frac{\partial (\rho w)}{\partial t} = - \nabla \cdot (\rho \mathbf{u} w) - \nabla p^\prime - \rho^\prime g + (\nabla \cdot \tau)_z + F^z,

where

.. math::

  p = \overline{p}(z) + p^\prime

and

.. math::

  \rho = \overline{\rho}(z) + \rho^\prime

and

.. math::

  \frac{d \overline{p}}{d z} = - \overline{\rho} g

with velocity :math:`\mathbf{u} = (u,v,w)` and gravity :math:`\mathbf{g} = (0,0,-g)`.

The relationship between potential temperature and temperature is given by

.. math::

  \theta = T (\frac{p_0}{p})^{R_d / c_p}

and we use the following equation of state:

.. math::

  p = \rho R_d T;

which can also be written in terms of :math:`\theta` as

.. math::

  p = (\rho R_d \theta / p_0^{R_d / c_p} )^\gamma

Here :math:`\rho, T, \theta`, and :math:`p` are the density, temperature, potential temperature and pressure, respectively;
these variables are all defined at cell centers.
:math:`A` is an advected quantity, i.e., a tracer, also defined at cell centers.
:math:`\mathbf{u}` and :math:`(\rho \mathbf{u})` are the velocity and momentum, respectively,
and are defined on faces.

:math:`R_d` and :math:`c_p` are the gas constant and specific heat capacity for dry air respectively,
and :math:`\gamma = c_p / (c_p - R_d)` .  :math:`p_0` is a reference value for pressure.

