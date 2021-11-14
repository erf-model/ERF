
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _Equations:



Compressible Navier-Stokes Equations
====================================

The following equations express conservation of mass, momentum, energy, and scalars in compressible fluid flow:

.. math::
  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}),

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + pI) +\rho \mathbf{g} + \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho \theta)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \theta) + \nabla \cdot (\alpha_{T}\ \nabla (\rho \theta)),

  \frac{\partial (\rho C)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} C) + \rho \nabla \cdot (\alpha_{C}\ \nabla C)

where :math:`\tau` is the stress tensor and :math:`\mathbf{F}` are the forcing terms described in :ref:`Forcings`.

The assumptions involved in deriving these equations from first principles are:

- Continuum behavior
- Ideal gas behavior (:math:`p = \rho R_d T`) with constant specific heats (:math:`c_p,c_v`)
- Constant mixture molecular weight (therefore constant :math:`R_d`)
- Viscous heating is negligible
- Newtonian fluid
- No chemical reactions, second order diffusive processes or radiative heat transfer
- Constant transport coefficients:  :math:`\mu, (\rho \alpha_T), (\rho \alpha_C)`. This is a good approximation for flows of
  interest because all are independent of density (or pressure), and only weakly dependent on temperature (:math:`T^{1/2}`)

.. note:: The diffusion coefficients, :math:`\alpha_{T}` and :math:`\alpha_{C}` for energy and scalar equations respectively,
   are in general variable for compressible flow. However, for low Mach number atmospheric flows they are assumed to be constant.

For low Mach number atmospheric flows the energy and scalar equations reduce to:

.. math::
  \frac{\partial (\rho \theta)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \theta) + \alpha_{T}\ \nabla^2 (\rho \theta)

  \frac{\partial (\rho C)}{\partial t}      &=& - \nabla \cdot (\rho \mathbf{u} C)      + \alpha_{C}\ \nabla^2 (\rho C)

DNS and LES
------------

Strain-rate tensor
~~~~~~~~~~~~~~~~~~
:math:`\sigma_{ij}`, is expressed as:

.. math::
   \sigma_{ij} = S_{ij} -D_{ij},

where 

.. math::
   S_{ij} = \frac{1}{2} \left(  \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}   \right)

.. math::
   D_{ij} = \frac{1}{3}  S_{kk} \delta_{ij} = \frac{1}{3} (\nabla \cdot \mathbf{u}) \delta_{ij}

DNS
~~~
When running in DNS mode, it is assumed that the mesh is fine enough to resolve the Kolmogorov scales of turbulence.
The dynamic viscosity (molecular), :math:`\mu`, is assumed to be constant for low Mach number atmospheric flows.

The stress tensor is thus given by 

.. math::
   \tau_{ij} = 2\mu \sigma_{ij}

.. note:: The bulk viscosity contribution to :math:`\tau_{ij}`, i.e., :math:`\kappa S_{kk} \delta_{ij}` has been ignored
   in the current implementation.

LES
~~~
When running in LES mode, it is assumed that the Kolmogorov scales are not resolved.  
LES models attempt to account for the effect of unresolved scales by replacing 
a constant :math:`\mu` by a turbulent viscosity, :math:`\mu_{t}`.

ERF offers two LES options: Smagorinsky and Deardorff models.

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

