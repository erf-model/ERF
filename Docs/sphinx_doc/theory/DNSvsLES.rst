
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
Therefore, in DNS mode the :ref:`previously described equations <GoverningEquations>` are solved directly with no additional models being required.
The transport coefficients correspond to the molecular transport coefficients,
with one of these assumptions:

* constant transport coefficients :math:`\mu`, :math:`\rho\alpha_C`, and :math:`\rho\alpha_T`, or

* constant :math:`\nu = \mu / \rho`, :math:`\alpha_T` and constant :math:`\alpha_C`,
      each of which is then multiplied by :math:`\rho` in the construction of the diffusive terms.

Note that scaling arguments from simple kinetic theory indicate that the quantities :math:`\mu`, :math:`\rho\alpha_C`, and :math:`\rho\alpha_T` are all independent of density and pressure but weakly dependent on temperature (:math:`T^{1/2}` scaling), justifying use of constant parameters over the range of conditions of interest.

LES
---
When running in Large Eddy Simulation (LES) mode, it is assumed that the Kolmogorov scales are not resolved. As a result, the numerical
discretization acts as a filter on the governing equations, resulting in the following set of filtered equations:

.. math::

  \frac{\partial \overline{\rho}}{\partial t} &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}}),

  \frac{\partial (\overline{\rho} \mathbf{\tilde{u}})}{\partial t} &= - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \mathbf{\tilde{u}}) - \nabla \overline{p} + \overline{\rho} \mathbf{g} - \nabla \cdot \overline{\tau} + \mathbf{\overline{F}} &- \nabla \cdot (\overline{\rho} \mathbf{\widetilde{u u}} - \overline{\rho}\mathbf{\tilde{u}\tilde{u}} ) ,

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

   \overline{\rho} {\widetilde{\mathbf{u} \theta}} - \overline{\rho}\mathbf{\tilde{u}}\tilde{\theta} &= -\frac{\mu_t}{Pr_t} \nabla \tilde{\theta}

   \overline{\rho} \widetilde{\mathbf{u} C} - \overline{\rho}\mathbf{\tilde{u}}\tilde{C} &= -\frac{\mu_t}{Sc_t} \nabla \tilde{C}

   \overline{\rho} \mathbf{\widetilde{u u}} - \overline{\rho}\mathbf{\tilde{u}\tilde{u}} &= \tau^{sfs}

.. math::

   \tau^{sfs}_{ij} - \frac{\delta_{ij}}{3} \tau^{sfs}_{kk} = -2 \mu_t \tilde{\sigma}_{ij}

   \tau^{sfs}_{kk} = 2 \mu_t \frac{C_I}{C_s^2} (2 \tilde{S}_{ij} \tilde{S}_{ij})^{1/2}.

The model coefficients :math:`C_s, C_I, Pr_t, Sc_t` have nominal values of 0.16, 0.09, 0.7, and 0.7,
respectively (Martin et al., Theoret. Comput. Fluid Dynamics (2000)).

.. note:: Presently, we assume :math:`C_I =0`. This term is similar to the bulk viscosity term for molecular transport and
      should be added if the bulk viscosity term is added. It is believed to be small for low-Mach number flows, but there
      is some discussion in the literature about this topic. See Moin et al., "A dynamic subgrid-scale model for
      compressible turbulence and scalar transport", PoF (1991); Martin et al., Subgrid-scale models for compressible
      large-eddy simulations", Theoret. Comput. Fluid Dynamics (2000).

When substituted back into the filtered equations, the gradient transport LES models take exactly the same form as the
molecular transport terms, but with the
constant molecular transport coefficients replaced by turbulent equivalents (e.g. :math:`\mu` becomes the turbulent viscosity,
:math:`\mu_{t}`). Therefore, when the code is run in LES mode, the :ref:`equation set<GoverningEquations>` remains the same,
but all variables are interpereted as the appropriate filtered version and the turbulent transport coefficients are used.
Molecular transport is omitted by default in the present LES implementation because the molecular
transport coefficients are insignificant compared to turbulent transport for most LES grids. However, for fine LES grids,
molecular transport and LES models may both be activated, with the effective transport coefficient being the sum of the molecular
and turbulent coefficients.

It should also be noted that filtering affects the computation of pressure from density and potential temperature, but the nonlinearity
in the equation of state is weak for :math:`\gamma = 1.4`, so the subfilter contribution is neglected:

.. math::
   \overline{p} = \overline{ \left( \frac{\rho R_d \theta}{p_0^{R_d / c_p}} \right)^\gamma} \approx \left( \frac{\overline{\rho} R_d \tilde{\theta}}{p_0^{R_d / c_p}} \right)^\gamma.

ERF offers two LES options: Smagorinsky and Deardorff models, which differ in how :math:`\mu_{t}` is computed.

.. _SmagorinskyModel:

Smagorinsky Model
~~~~~~~~~~~~~~~~~~

.. math::
   \mu_{t} = (C_s \Delta)^2 (\sqrt{2 \tilde{S} \tilde{S}}) \overline{\rho}

:math:`C_s` is the Smagorinsky constant and :math:`\Delta` is the cube root of cell volume, the representative mesh spacing.

.. math::
   \tau_{ij} = -2\mu_{t} \tilde{\sigma}_{ij} = -K \tilde{\sigma}_{ij}

where :math:`K = 2\mu_{t}`

In the Smagorinsky model, modeling of :math:`\mu_{t}` does not account for the turbulent kinetic energy (TKE) corresponding to
unresolved scales and no extra equation for TKE is solved.

Enhanced Smagorinsky Model with Moist Stratification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The standard Smagorinsky model uses a fixed relationship between mixing length and grid scale (:math:`\ell = \Delta`).
For atmospheric applications involving moisture, this can be inappropriate in stably stratified regions where buoyancy effects
significantly suppress turbulent mixing. ERF includes an enhanced Smagorinsky formulation that accounts for moist atmospheric
stratification effects while maintaining computational efficiency.

The enhanced model modifies the mixing length based on local atmospheric stability:

.. math::
   \ell_m = \begin{cases}
   \Delta & \text{if } N^2_{eff} \leq 0 \\
   \min\left(\Delta, \sqrt{\frac{u_c^2}{N^2_{eff}}}\right) & \text{if } N^2_{eff} > 0
   \end{cases}

where :math:`u_c = |\tilde{S}| \Delta` is a characteristic velocity scale based on the local strain rate magnitude,
and :math:`N^2_{eff}` is the effective moist Brunt-Väisälä frequency.

**Moist Stratification Parameter**

The moist stratification parameter accounts for three key atmospheric effects:

1. **Virtual potential temperature effects**: Water vapor decreases air density (enhancing buoyancy) while condensate increases density (reducing buoyancy):

.. math::
   \theta_v = \theta \left(1 + 0.61 q_v - q_c - q_i\right)

where :math:`q_v`, :math:`q_c`, and :math:`q_i` are mixing ratios of water vapor, liquid water, and ice, respectively.

2. **Moist Brunt-Väisälä frequency**: The stratification parameter uses virtual potential temperature instead of dry potential temperature:

.. math::
   N^2_m = \frac{g}{\theta_0} \frac{\partial \theta_v}{\partial z}

3. **Conditional instability correction**: In saturated regions with condensate present, the effective stratification is reduced due to latent heat release:

.. math::
   N^2_{eff} = N^2_m \cdot \frac{1 + \frac{L_{eff} q_{s,eff}}{R_d T}}{1 + \frac{L_{eff}^2 q_{s,eff}}{c_p R_v T^2}}

where :math:`L_{eff}` and :math:`q_{s,eff}` are phase-weighted effective latent heat and saturation mixing ratio, accounting for the presence of both liquid and ice phases:

.. math::
   L_{eff} &= f_l L_v + f_i (L_v + L_i) \\
   q_{s,eff} &= f_l q_{s,w} + f_i q_{s,i} \\
   f_l &= \frac{q_c}{q_c + q_i + \epsilon}, \quad f_i = 1 - f_l

The enhanced turbulent viscosity becomes:

.. math::
   \mu_{t} = (C_s \ell_m)^2 (\sqrt{2 \tilde{S} \tilde{S}}) \overline{\rho}

with mixing length bounds applied to prevent numerical issues:

.. math::
   \ell_m \geq 0.001 \Delta

**Physical Interpretation**

This formulation creates an effective Richardson number criterion where turbulent length scales respond to atmospheric stability:

.. math::
   Ri_{eff} = \frac{N^2_{eff}}{|\tilde{S}|^2}

When :math:`Ri_{eff} > 1`, stratification dominates and mixing is suppressed.
When :math:`0 < Ri_{eff} < 1`, shear dominates and mixing approaches standard Smagorinsky behavior.
When :math:`Ri_{eff} \leq 0`, stratification is unstable and the standard Smagorinsky behavior is exactly recovered, regardless of moisture.

**Implementation Control**

The moist Richardson number limiter can be enabled via the input parameter:

.. code-block:: none

    erf.use_moist_Ri_correction = true   # Apply moist Richardson number limiter
    erf.use_moist_Ri_correction = false  # Standard Smagorinsky

When disabled, the model reverts to the standard Smagorinsky formulation with :math:`\ell_m = \Delta`.

Deardorff Model
~~~~~~~~~~~~~~~

Unlike the Smagorinsky model, the Deardorff model accounts for the contribution of TKE in modeling :math:`\mu_{t}` and a prognostic equation
for TKE is solved.  The turbulent viscosity is computed as:

.. math::
   \mu_t = \overline{\rho} C_k \ell (e^{sfs})^{1/2},

where the mixing length :math:`\ell = \Delta` for unstable (or neutral) stratification, otherwise the mixing length is reduced as per Deardorff 1980:

.. math::
   \ell = \frac{0.76 (e^{sfs})^{1/2}}{\left(\frac{g}{\theta_0}\frac{\partial\theta}{\partial z}\right)^{1/2}}.

The potential temperature gradient in the denominator dictates the stratification, which is scaled by a reference virtual potential temperature :math:`\theta_0`.
The mixing length is set to a maximum of :math:`\Delta` to prevent unbounded behavior under weakly stable conditions.

Then the equation solved to determine :math:`e^{sfs}`, the subfilter contribution to TKE, is:

.. math::
   \frac{\partial \overline{\rho} e^{sfs}}{\partial t} = - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{e}^{sfs})
                                                         - \tau_{ij} \frac{\partial \tilde{u}_i}{\partial x_j}
                                                         + \frac{g}{\theta_0} \tau_{\theta w}
                                                         + \nabla \cdot \left( \frac{\mu_t}{\sigma_k} \nabla e^{sfs}  \right)
                                                         - \overline{\rho} C_\epsilon \frac{(e^{sfs})^{3/2}}{\overline{\Delta}}.

where :math:`\sigma_k` is a constant model coefficient representing the ratio of turbulent viscosity to turbulent diffusivity
of TKE that should be order unity (e.g., Moeng 1984 uses TKE diffusivity of :math:`2 \mu_t`), we have used the downgradient diffusion assumption

.. math::
   \frac{\partial\left\langle \left( u_{n}^{'}\rho e + u_{n}^{'}p^{'} \right) \right\rangle}{\partial x_{n}} =
           -\nabla \cdot \left( \frac{\mu_t}{\sigma_k} \nabla e^{sfs}  \right),

the eddy diffusivity of heat is

.. math::
   K_H = \left(1 + \frac{2\ell}{\Delta}\right) \mu_t,

and the SFS heat flux is

.. math::
   \tau_{\theta i} = -K_H \frac{\partial\theta}{\partial x_i}.

The RHS terms of the TKE transport equation correspond to advection, shear production, buoyant production, diffusion, and dissipation.
