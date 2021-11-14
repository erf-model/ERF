
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

  \frac{\partial (\rho C)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} C) + \nabla \cdot (\alpha_{C}\ \nabla (\rho C))

Here :math:`\tau_{ij}` is the stress tensor and is related to shear-rate tensor, :math:`\sigma_{ij}`,  by:

.. math::
   \tau_{ij} = 2\mu\sigma_{ij}

The shear-rate tensor, :math:`\sigma_{ij}`, is further expressed as:

.. math::
   \sigma_{ij} = S_{ij} -D_{ij},

where :math:`S_{ij}` is the strain-rate tensor and :math:`D_{ij}` is the expansion-rate tensor.

.. math::
   S_{ij} = \frac{1}{2} \left(  \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}   \right)

.. math::
   D_{ij} = \frac{1}{3}  S_{kk} \delta_{ij} = \frac{1}{3} (\nabla \cdot \mathbf{u}) \delta_{ij}

.. note:: For an incompressible flow :math:`\nabla \cdot \mathbf{u} = 0`,
   where as for compressible flow :math:`\nabla \cdot \mathbf{u} \neq 0` in general,
   thus leading to a non-negligible :math:`D_{ij}`.

This leads to:

.. math::
   \tau_{ij} = 2\mu \left( S_{ij} - \frac{1}{3} S_{kk} \delta_{ij} \right), \hspace{24pt}

.. note:: The bulk viscosity contribution to :math:`\tau_{ij}`, i.e., :math:`\kappa S_{kk} \delta_{ij}` has been ignored
   in the current implementation.

.. note:: The diffusion coefficients, :math:`\alpha_{T}` and :math:`\alpha_{C}` for energy and scalar equations respectively,
   are in general variable for compressible flow. However, for low Mach number atmospheric flows they are assumed to be constant.

For low Mach number atmospheric flows the energy and scalar equations reduce to:

.. math::
  \frac{\partial (\rho \theta)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \theta) + \alpha_{T}\ \nabla^2 (\rho \theta)

  \frac{\partial (\rho C)}{\partial t}      &=& - \nabla \cdot (\rho \mathbf{u} C)      + \alpha_{C}\ \nabla^2 (\rho C)

DNS and LES
------------

DNS
~~~
When run in DNS mode, it is assumed that the mesh is fine enough to resolve Kolmogorov scales of turbulence.
The viscosity, :math:`\mu`, used in computing the stress term, :math:`\tau_{ij}`, is the dynamic viscosity (molecular).
:math:`\mu` is assumed to be constant for low Mach number atmospheric flows and is provided by the user.

LES
~~~
When using the LES model, it is assumed that mesh is coarser than those used for DNS and Kolmogorov scales are not resolved.
The scales that are resolved are represented by the filtered equations and the unresolved scales are modeled.
The LES models essentially attempt to account for the effect of unresolved scales on the viscosity where
:math:`\mu` is replaced by a modeled turbulent viscosity, :math:`\mu_{t}`.

There are several approaches to model :math:`\mu_{t}`. Smagorinsky and Deardorff models are commonly used.

Smagorinsky Model
~~~~~~~~~~~~~~~~~~
.. math::
   \mu_{t} = (C_s \Delta)^2 (\sqrt{2 S S}) \rho
:math:`C_s` is the Smagorinsky constant and :math:`\Delta` is the cube root of cell volume, the representative mesh spacing.

.. math::
   \tau_{ij} = 2\mu_{t} \sigma_{ij} = K \sigma_{ij}

where :math:`K = 2\mu_{t}`

In Smagorisnky model, modeling of :math:`\mu_{t}` doesn't account for the turbulent kinetic energy (TKE) corresponding to
unresolved scales and no extra equation for TKE is solved.

Deardorff Model
~~~~~~~~~~~~~~~~~~
Unlike Smagorinsky model, Deardorff model accounts for the contribution of TKE in modeling :math:`\mu_{t}` and equation
for TKE is solved, i.e., TKE is prognosed.

Forcings
------------
:math:`\mathbf{F}` are the forcing terms described in :ref:`Forcings`.


Equations in Perturbation Form
-------------------------------
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

Diagnostic Relationships
-------------------------

The relationship between potential temperature and temperature is given by

.. math::

  \theta = T \left( \frac{p_0}{p} \right)^{R_d / c_p}

and we use the following equation of state:

.. math::

  p = \rho R_d T;

which can also be written in terms of :math:`\theta` as

.. math::

  p = \left( \frac{\rho R_d \theta}{p_0^{R_d / c_p}} \right)^\gamma

Here :math:`\rho, T, \theta`, and :math:`p` are the density, temperature, potential temperature and pressure, respectively;
these variables are all defined at cell centers.
:math:`A` is an advected quantity, i.e., a tracer, also defined at cell centers.
:math:`\mathbf{u}` and :math:`(\rho \mathbf{u})` are the velocity and momentum, respectively,
and are defined on faces.

:math:`R_d` and :math:`c_p` are the gas constant and specific heat capacity for dry air respectively,
and :math:`\gamma = c_p / (c_p - R_d)` .  :math:`p_0` is a reference value for pressure.

