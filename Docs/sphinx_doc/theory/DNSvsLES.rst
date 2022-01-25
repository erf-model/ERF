
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _DNSvsLES:

Simulation Modes: DNS and LES
=============================

DNS
---

When running in Direct Numerical Simulation (DNS) mode, it is assumed that the mesh is fine enough to resolve the Kolmogorov scales of turbulence.
Therefore, in DNS mode the equations above are solved directly with no additional models being required.

LES
---
When running in Large Eddy Simulation (LES) mode, it is assumed that the Kolmogorov scales are not resolved. As a result, the numerical
discretization acts as a filter on the governing equations, resulting in the following set of filtered equations:

.. math::

  \frac{\partial \overline{\rho}}{\partial t} &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}}),

  \frac{\partial (\overline{\rho} \mathbf{\tilde{u}})}{\partial t} &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \mathbf{\tilde{u}}) - \nabla \overline{p} + \overline{\rho} \mathbf{g} + \nabla \cdot \overline{\tau} + \mathbf{\overline{F}} &- \nabla \cdot (\overline{\rho} \mathbf{\widetilde{u u}} - \overline{\rho}\mathbf{\tilde{u}\tilde{u}} ) ,

  \frac{\partial (\overline{\rho} \tilde{\theta})}{\partial t} &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{\theta}) + \nabla \cdot \left( \overline{\rho} \alpha_{T} \nabla \tilde{\theta} \right)  &- \nabla \cdot (\overline{\rho} {\widetilde{\mathbf{u} \theta}} - \overline{\rho}\mathbf{\tilde{u}}\tilde{\theta} ) ,

  \frac{\partial (\overline{\rho} \tilde{C})}{\partial t}      &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{C})      + \nabla \cdot \left( \overline{\rho} \alpha_{C} \nabla \tilde{C} \right)  &- \nabla \cdot (\overline{\rho} \widetilde{\mathbf{u} C} - \overline{\rho}\mathbf{\tilde{u}}\tilde{C} ) ,

where overbars indicate filtering and tildes indicate density-weighted (Favre) filtering
(e.g., :math:`\tilde{\theta} = \overline{\rho \theta} / \overline{\rho}`).
When the code is run in LES mode, all variables correspond to their appropriate filtered version.

In the above equations, the final term in each of the momentum, potential temperature, and scalar equations is unclosed
due to containing a filtered nonlinear function of the state quantities. These terms represent the effect of turbulent transport at unresolved scales.
LES models attempt to account for these terms by
invoking a gradient transport hypothesis, which assumes that turbulent transport acts similarly to molecular transport
in that quantities are transported down their resolved gradients:

.. math::

   \overline{\rho} {\widetilde{\mathbf{u} \theta}} - \overline{\rho}\mathbf{\tilde{u}}\tilde{\theta} &= -frac{\mu_t}{Pr_t} \nabla \tilde{\theta}

   \overline{\rho} \widetilde{\mathbf{u} C} - \overline{\rho}\mathbf{\tilde{u}}\tilde{C} &= -frac{\mu_t}{Sc_t} \nabla \tilde{C}

   \overline{\rho} \mathbf{\widetilde{u u}} - \overline{\rho}\mathbf{\tilde{u}\tilde{u}} &= -\tau^{sfs}

.. math::

   \tau^{sfs}_{ij} - \frac{\delta_{ij}}{3} \tau^{sfs}_{kk} = 2 \mu_t \tilde{\sigma}_{ij}

   \tau^{sfs}_{kk} = 2 \mu_t \frac{C_I}{C_s^2} (2 \tilde{S}_{ij} \tilde{S}_{ij})^{1/2}.

The model coefficients :math:`C_s, C_I, Pr_t, Sc_t` have nominal values of 0.16, 0.09, 0.7, and 0.7,
respectively (Martin et al., Theoret. Comput. Fluid Dynamics (2000)).
Note that the gradient transport LES models take exactly the same form as the molecular transport terms, but with the
constant molecular transport coefficients replaced by turbulent equivalents (e.g. :math:`\mu` becomes the turbulent viscosity,
:math:`\mu_{t}`). Molecular transport is omitted by default in the present implementation because the molecular
transport coefficients are insignificant compared to turbulent transport for most LES grids. However, for fine LES grids, molecular transport and LES models may both be activated.

.. note:: Presently, we assume :math:`C_I =0`. This term is similar to the bulk viscosity term for molecular transport and
      should be added if the bulk viscosity term is added. It is believed to be small for low-Mach number flows, but there
      is some discussion in the literature about this topic. See Moin et al., "A dynamic subgrid-scale model for
      compressible turbulence and scalar transport", PoF (1991); Martin et al., Subgrid-scale models for compressible
      large-eddy simulations", Theoret. Comput. Fluid Dynamics (2000).

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

Then the equation solved to determine :math:`e^{sfs}`, the subfilter contribution to TKE, is:

.. math::

   \frac{\partial \overline{\rho} e^{sfs}}{\partial t} = - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{e}^{sfs})
                                                         + K_H \rho (\frac{g}{\theta} \frac{\partial\theta}{\partial z})
                                                         + \nabla \cdot \left( \frac{\mu_t}{\sigma_k} \nabla e^{sfs}  \right)
                                                         - \tau_{ij} \frac{\partial \tilde{u}_i}{\partial x_j}
                                                         - \overline{\rho} C_\epsilon \frac{(e^{sfs})^{3/2}}{\overline{\Delta}}.

where :math:`\sigma_k` is a constant model coefficient representing the ratio of turbulent viscosity to turbulent diffusivity
of TKE that should be order unity, and we have used

.. math::

   \frac{\partial\left\langle \left( u_{n}^{'}\rho e + u_{n}^{'}p^{'} \right) \right\rangle}{\partial x_{n}} =
           -\nabla \cdot \left( \frac{\mu_t}{\sigma_k} \nabla e^{sfs}  \right)
