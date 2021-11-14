
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _Equations:



Compressible Navier-Stokes Equations
====================================

The following partial differential equations (prognostic equations) express conservation of mass, momentum, energy, and scalars in compressible fluid flow:

.. math::
  \frac{\partial \rho}{\partial t} &= - \nabla \cdot (\rho \mathbf{u}),

  \frac{\partial (\rho \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + pI) +\rho \mathbf{g} + \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho \theta)}{\partial t} &= - \nabla \cdot (\rho \mathbf{u} \theta) + \nabla \cdot (\alpha_{T}\ \nabla (\rho \theta)),

  \frac{\partial (\rho C)}{\partial t} &= - \nabla \cdot (\rho \mathbf{u} C) +\nabla \cdot (\rho \alpha_{C}\ \nabla C)

where :math:`\tau` is the viscous stress tensor and :math:`\mathbf{F}` are the forcing terms described in :ref:`Forcings`, and the potential temperature :math:`\theta` is defined

.. math::
   
  \theta = T \left( \frac{p_0}{p} \right)^{R_d / c_p}.

The assumptions involved in deriving these equations from first principles are:

- Continuum behavior
- Ideal gas behavior (:math:`p = \rho R_d T`) with constant specific heats (:math:`c_p,c_v`)
- Constant mixture molecular weight (therefore constant :math:`R_d`)
- Viscous heating is negligible
- No chemical reactions, second order diffusive processes or radiative heat transfer

We further assume constant transport coefficients :math:`\mu`, :math:`(\rho \alpha_C)`, and :math:`(\alpha_T)`.
This is a good approximation for flows of interest because all are independent of density (or pressure),
and only weakly dependent on temperature (:math:`T^{1/2}` scaling based on kinetic theory).  Thus for
low Mach number atmospheric flows the energy and scalar equations reduce to:

.. math::
  \frac{\partial (\rho \theta)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \theta) + \nabla \cdot (\alpha_{T}\ \nabla (\rho \theta)),


  \frac{\partial (\rho C)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} C) + \rho \nabla \cdot (\alpha_{C}\ \nabla C)

Similarly, the viscous stress tensor can be expressed:

.. math:: 
   \tau_{ij} = 2\mu \sigma_{ij}

where the anisotropic strain rate tensor :math:`\sigma_{ij} = S_{ij} -D_{ij}`,

.. math::
   S_{ij} = \frac{1}{2} \left(  \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}   \right), \hspace{24pt}
   D_{ij} = \frac{1}{3}  S_{kk} \delta_{ij} = \frac{1}{3} (\nabla \cdot \mathbf{u}) \delta_{ij},

and the  bulk viscosity contribution to :math:`\tau_{ij}`, i.e., :math:`\kappa S_{kk} \delta_{ij}` has been ignored
in the current implementation.

Diagnostic Relationships
------------------------

In order to close the above prognostic equations, a relationship between the pressure and the reesnported state variables
must be specified. This is obtained by re-expressing the ideal gas equation of state in terms of :math:`\theta`:

.. math::
   p = \left( \frac{\rho R_d \theta}{p_0^{R_d / c_p}} \right)^\gamma

Nomenclature
------------
Here :math:`\rho, T, \theta`, and :math:`p` are the density, temperature, potential temperature and pressure, respectively;
these variables are all defined at cell centers.
:math:`C` is an advected quantity, i.e., a tracer, also defined at cell centers.
:math:`\mathbf{u}` and :math:`(\rho \mathbf{u})` are the velocity and momentum, respectively,
and are defined on faces.

:math:`R_d` and :math:`c_p` are the gas constant and specific heat capacity for dry air respectively,
and :math:`\gamma = c_p / (c_p - R_d)` .  :math:`p_0` is a reference value for pressure.


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

Simulation Modes: DNS and LES
=============================
  
DNS
---

When running in Direct Numerical Simulation (DNS) mode, it is assumed that the mesh is fine enough to resolve the Kolmogorov scales of turbulence.
Therefore, in DNS mode the equations above are solved directly with no additional models being required.

LES
---
When running in Large Eddy SImulation (LES) mode, it is assumed that the Kolmogorov scales are not resolved. As a result, the numerical
discretization acts as a filter on the governing equations, resulting in the following set of filtered equations:

.. math::

  \frac{\partial \overline{\rho}}{\partial t} &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}}),

  \frac{\partial (\overline{\rho} \mathbf{\tilde{u}})}{\partial t} &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \mathbf{\tilde{u}} + \overline{p}I) +\overline{\rho} \mathbf{g} + \nabla \cdot \overline{\tau} + \mathbf{\overline{F}} &- \nabla \cdot (\overline{\rho} \mathbf{\widetilde{u u}} - \overline{\rho}\mathbf{\tilde{u}\tilde{u}} ) ,

  \frac{\partial (\overline{\rho} \tilde{\theta})}{\partial t} &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{\theta}) + \overline{\rho} \alpha_{T} \nabla^2 \tilde{\theta}  &- \nabla \cdot (\overline{\rho} {\widetilde{\mathbf{u} \theta}} - \overline{\rho}\mathbf{\tilde{u}}\tilde{\theta} ) ,

  \frac{\partial (\overline{\rho} \tilde{C})}{\partial t}      &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{C})      + \overline{\rho} \alpha_{C} \nabla^2 \tilde{C}  &- \nabla \cdot (\overline{\rho} \widetilde{\mathbf{u} C} - \overline{\rho}\mathbf{\tilde{u}}\tilde{C} ) ,

where overbars indicate filtering and tildes indicate density-weighted (Favre) filtering
(e.g., :math:`\tilde{\theta} = \overline{\rho \theta} / \overline{\rho}`).
When the code is run in LES mode, all variables correspond to their appropriate filtered version. 

In the above equations, the final term in each of the momentum, potential temperature, and scalar equations is unclosed
due to containing a filtered nonlinear function of the state quantities. These terms represent the effect of turbulent transport at unresolved scales. 
LES models attempt to account for these terms by
invoking a gradient transport hypothesis, which assumes that turbulent transport acts similarly to molecular transport
in that quantities are transported down their resolved gradients:

.. math::

   \overline{\rho} {\widetilde{\mathbf{u} \theta}} - \overline{\rho}\mathbf{\tilde{u}}\tilde{\theta} &= \frac{\mu_t}{Pr_t} \nabla \tilde{\theta}

   \overline{\rho} \widetilde{\mathbf{u} C} - \overline{\rho}\mathbf{\tilde{u}}\tilde{C} &= \frac{\mu_t}{Sc_t} \nabla \tilde{C}

   \overline{\rho} \mathbf{\widetilde{u u}} - \overline{\rho}\mathbf{\tilde{u}\tilde{u}}  &=  \tau^{sfs}

.. math::

   \tau^{sfs}_{ij} - \frac{\delta_{ij}}{3} \tau^{sfs}_{kk} = 2 \mu_t \tilde{\sigma}_{ij}

   \tau^{sfs}_{kk} = 2 \mu_t \frac{C_I}{C_s^2} (2 \tilde{S}_{ij} \tilde{S}_{ij})^{1/2}.
   
The model coefficients :math:`C_s, C_I, Pr_t, Sc_t` have nominal values of 0.16, 0.09, 0.7, amd 0.7,
respectively (Martin et al., Theoret. Comput. Fluid Dynamics (2000)).
Note that the gradient transport LES models take exactly the same form as the molecular transport terms, but with the constant
constant molecular transport coefficients replaced by turbulent requivalents (e.g. :math:`\mu` becomes the turbulent viscosity,
:math:`\mu_{t}`). Molecular transport is omitted in the present implementation because the molecular
transport coefficients are insignificant compared to turbulent transport for most LES grids.

.. note:: The omission of molecular transport in LES mode will need to be revisisted if resolutions close to DNS become of interest.
      Presently, we also assume :math:`C_I =0`. This term is similar to the bulk viscosity term for molecular transport and
      should be added if the bulk viscosity term is added. It is beieved to be small for low-Mach number flows, but there
      is some discussion in the literature about this topic. See Moin et al., "A dynamic subgrid-scale model for
      compressible turbulence and scalar transport", PoF (1991); Martin et al., Subgrid-scale models for compressible
      large-eddy simulations", Theoret. Comput. Fluid Dynamics (2000).

.. note:: LES models for potential temperature and scalars have not yet been implemented, and the molecular transport coefficients
      are retained in these equations. This will be updated soon.

It should also be noted that filtering affects the computation of pressure from density and potential temperature, but the nonlinearity
in the equation of state is weak for :math:`\gamma = 1.4`, so the subfilter contribution is neglected:

.. math::
   \overline{p} = \overline{ \left( \frac{\rho R_d \theta}{p_0^{R_d / c_p}} \right)^\gamma} \approx \left( \frac{\overline{\rho} R_d \tilde{\theta}}{p_0^{R_d / c_p}} \right)^\gamma.

ERF offers two LES options: Smagorinsky and Deardorff models, which differ in how :math:`\mu_{t}` is computed.

Smagorinsky Model
~~~~~~~~~~~~~~~~~~
.. math::
   \mu_{t} = (C_s \Delta)^2 (\sqrt{2 \tilde{S} \tilde{S}}) \overline{\rho}
:math:`C_s` is the Smagorinsky constant and :math:`\Delta` is the cube root of cell volume, the representative mesh spacing.

.. math::
   \tau_{ij} = 2\mu_{t} \tilde{\sigma}_{ij} = K \tilde{\sigma}_{ij}

where :math:`K = 2\mu_{t}`

In the Smagorinsky model, modeling of :math:`\mu_{t}` does not account for the turbulent kinetic energy (TKE) corresponding to
unresolved scales and no extra equation for TKE is solved.

Deardorff Model
~~~~~~~~~~~~~~~
Unlike the Smagorinsky model, the Deardorff model accounts for the contribution of TKE in modeling :math:`\mu_{t}` and a prognostic equation
for TKE is solved.  The turbulent viscosity is computed as:

.. math::

   \mu_t = C_k \overline{\rho} \Delta (k^{sfs})^{1/2}.

The equation solved to determine :math:`k^{sfs}`, the subfilter contribution to TKE, is:
   
.. math::

   \frac{\partial \overline{\rho} k^{sfs}}{\partial t} = - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{k}^{sfs}) + \nabla \cdot \left( \frac{\mu_t}{\sigma_k} \nabla k ^{sfs}  \right) + ( \overline{\rho} \widetilde{\mathbf{uu}} - \overline{\rho} \tilde{\mathbf{u}} \tilde{\mathbf{u}})\nabla \cdot \mathbf{\tilde{u}} - \overline{\rho} C_\epsilon \frac{(k^{sfs})^{3/2}}{\overline{\Delta}}.

