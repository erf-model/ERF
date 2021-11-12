
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _Equations:



Compressible Navier-Stokes Equations
====================================


The following set of equations express conservation of mass, momentum, energy, and scalars in compressible fluid flow:

.. math::

  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}),

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + pI) +\rho \mathbf{g} + \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho T)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} T) + \rho \alpha_{T} \nabla^2 T + \frac{1}{c_p}\frac{Dp}{Dt},

  \frac{\partial (\rho C)}{\partial t}      &=& - \nabla \cdot (\rho \mathbf{u} C)      + \rho \alpha_{C} \nabla^2 C,

.. note:: When I look in the source code, it seems like we are solving diffusion as :math:`\alpha_{C} \nabla^2 \rho C` (i.e., similar to what is in the documentation elsewhere: :ref:`AcousticSubstep`), not what is written above. Is there something I'm missing? I believe what is written here (and in the discretization section) is what we should be doing. Regarding the transport coefficients, the fundamental transport coefficients :math:`\mu, (\rho \alpha_T), (\rho \alpha_C)` are independent of density and pressure and scale with the square root of temperature, so in cases of interest should be roughly constant - but note this implies inverse dependence of :math:`\alpha_i` on density. In any event, constant (:math:`\rho \alpha_C`), etc. is consistent with the above equations. 

where

.. math::

   \tau_{ij} = 2\mu \left( S_{ij} - \frac{2}{3} S_{kk} \delta_{ij} \right), \hspace{24pt}  S_{ij} = \frac{1}{2} \left(  \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}   \right)

.. note:: ERF is missing the subtraction of the isotropic part of the strain rate. For full generality, we would want to add the bulk viscosity contribution to :math:`\tau_{ij}` as well (:math:`\kappa S_{kk} \delta_{ij}`), but it is likely negligible for our cases, so this is not essential.
   
.. \tau = \mu \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}  - \frac{2}{3}\frac{\partial u_k}{\partial x_k} \delta_{ij} \right)
   
The assumptions involved in deriving these equations from first principles are:

- Continuum behavior
- Ideal gas behavior (:math:`p = \rho R_d T`) with constant specific heats (:math:`c_p,c_v`)
- Constant mixture molecular weight (therefore constant :math:`R_d`)
- Viscous heating is negligible
- Newtonian fluid
- No chemical reaction
- No second order diffusive processes
- No radiative heat transfer
- Stress due to the bulk viscosity is negligible
- Constant transport coefficients:  :math:`\mu, (\rho \alpha_T), (\rho \alpha_C)`. This is a good approximation for flows of
  interest because all are independent of density (or pressure), and only weakly dependent on temperature (:math:`T^{1/2}`)
- The body forces :math:`\mathbf{F}` are the forcing terms described in :ref:`Forcings`

.. note:: This list is meant to be exhaustive. What is missing?  

DNS
---

In DNS mode, ERF advances the above equations with the exception that the temperature :math:`T` is replaced by the potential temperature :math:`\theta`:

.. note:: Someone who understands atmospheric science better than I do should add a sentence here explaining why we use :math:`\theta`.
	  
.. math::

   \theta = T \left( \frac{p_0}{p} \right)^{R_d / c_p},

with the porgnostic (transport) equation

.. math::

  \frac{\partial (\rho \theta)}{\partial t} = - \nabla \cdot (\rho \mathbf{u} \theta) + \rho \alpha_{T} \nabla^2 \theta,

.. note:: When I tried to re-derive this equation from the temperature equation above, I arrived at :math:`\rho \alpha_{T} \left(\frac{p_0}{p}\right)^{R_d/c_p} \nabla^2 \theta` for the conduction term, and could not find a way to further simplify to the above form. How is that term obtained?  
  
which results in a modified equation of state

.. math::

     p = \left( \frac{\rho R_d \theta}{p_0^{R_d / c_p}} \right)^\gamma.


LES
---

     
Filtering the above equations results in the governing equations to be used in LES:

.. math::
   
  \frac{\partial \overline{\rho}}{\partial t} &=& - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}}),

  \frac{\partial (\overline{\rho} \mathbf{\tilde{u}})}{\partial t} &=& - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \mathbf{\tilde{u}} + \overline{p}I) +\overline{\rho} \mathbf{g} + \nabla \cdot \overline{\tau} + \mathbf{\overline{F}} &- \nabla \cdot (\overline{\rho} \mathbf{\widetilde{u u}} - \overline{\rho}\mathbf{\tilde{u}\tilde{u}} ) ,

  \frac{\partial (\overline{\rho} \tilde{\theta})}{\partial t} &=& - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{\theta}) + \overline{\rho} \alpha_{T} \nabla^2 \tilde{\theta}  &- \nabla \cdot (\overline{\rho} {\widetilde{\mathbf{u} \theta}} - \overline{\rho}\mathbf{\tilde{u}}\tilde{\theta} ) ,

  \frac{\partial (\overline{\rho} \tilde{C})}{\partial t}      &=& - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{C})      + \overline{\rho} \alpha_{C} \nabla^2 \tilde{C}  &- \nabla \cdot (\overline{\rho} \widetilde{\mathbf{u} C} - \overline{\rho}\mathbf{\tilde{u}}\tilde{C} ) ,

.. note:: The molecular transport terms are shown being retained here, but not in the present implementation. While they are likely to be small reative to the turbulent transport terms, I don't see a downside to keeping them around, or at least have the option to do so.
	  
where 

.. math::

   \overline{\tau}_{ij} = 2\mu \left( \overline{S}_{ij} - \frac{2}{3} \overline{S}_{kk} \delta_{ij} \right), \hspace{24pt}  \overline{S}_{ij} = \frac{1}{2} \left(  \frac{\partial \overline{u}_i}{\partial x_j} + \frac{\partial \overline{u}_j}{\partial x_i}   \right)

.. math::

   \overline{p} = \overline{ \left( \frac{\rho R_d \theta}{p_0^{R_d / c_p}} \right)^\gamma} . 

Any term that is a nonlinear function of the Favre-filtered state variables is in principal unclosed and requires some form of model. Additionally, for completeness, the viscous term depends on non-density-weighted filtered velocity (:math:`\mathbf{\overline{u}}`), while the equation is solved for density-weighted-filtered velocity (:math:`\mathbf{\tilde{u}}`). However, because the importance of this term is small, and extreme density variation over single grid cells is unlikely, computing this term using the density-weighted filtered velocity should not result in significant error relative to other models.  Based on the assumptions for this system, the terms that require modeling include the turbulent transport terms (final term in each of the momentum, potential temperature, and scalar equations) and the equation of state. For :math:`\gamma=1.4`, the equation of state is only weakly nonlinear, and neglecting the nonlinearity is likely yield satisfactory results:

.. math::

   \overline{p} \approx  \left( \frac{\overline{\rho} R_d \tilde{\theta}}{p_0^{R_d / c_p}} \right)^\gamma.

The turbulent transport terms can be closed with eddy viscosity-based models:

.. math::

   \overline{\rho} {\widetilde{\mathbf{u} \theta}} - \overline{\rho}\mathbf{\tilde{u}}\tilde{\theta} &=& \tau^{sfs}_{\theta} &=& \overline{\rho} \frac{\nu_t}{Pr_t} \nabla \theta

   \overline{\rho} \widetilde{\mathbf{u} C} - \overline{\rho}\mathbf{\tilde{u}}\tilde{C} &=&  \tau^{sfs}_{C}  &=&  \overline{\rho} \frac{\nu_t}{Sc_t} \nabla C

   \overline{\rho} \mathbf{\widetilde{u u}} - \overline{\rho}\mathbf{\tilde{u}\tilde{u}}  &=& \tau^{sfs}_{ij}

.. math::
   
   \tau^{sfs}_{ij} - \frac{\delta_{ij}}{3} \tau^{sfs}_{kk} = 2\overline{\rho} \nu_t \left( \tilde{S}_{ij} - \frac{2}{3} \tilde{S}_{kk} \delta_{ij} \right)

   \tau^{sfs}_{kk} = 2 \overline{\rho} \nu_t \frac{C_I}{C_s^2} (2 S_{ij} S_{ij})^{1/2}

.. note::
   The last closure for the isotropic part of the subfilter stresses may be unnecessary, particularly for flows that
   have low Mach numbers, i.e. :math:`C_I = 0`. See discussion in Moin et al., "A dynamic subgrid-scale model for
   compressible turbulence and scalar transport", PoF (1991); Martin et al., Subgrid-scale models for compressible
   large-eddy simulations", Theoret. Comput. Fluid Dynamics (2000). Right now, we do not model the isotropic portion
   separately (:math:`\tau^{sfs}_{ij} = 2\overline{\rho} \nu_t \tilde{S}_{ij}`. However, the term is modelled in
   PeleC and some remnants of that seems to have made their way into ERF (:math:`C_I` exists as a parameter)

We consider two types of LES models, Smagorinsky and Deardorff, which differ in how :math:`\nu_t` is computed.

Smagorinsky
~~~~~~~~~~~
   
In Smagorinsky models, the turbulent viscosity is approximated with an algebraic/diagnostic equation as:

.. math:: \nu_t = (C_s \overline{\Delta})^2 (2 S_{ij} S_{ij})^{1/2}

The model coefficients :math:`C_s, C_I, Pr_t, Sc_t` have nominal values of 0.16, 0.09, 0.7, amd 0.7, respectively (Martin et al., Theoret. Comput. Fluid Dynamics (2000)), but may also be determined by a dynamic procedure (Moin et al., PoF (1991)). 

.. note:: Do we want to implement a dynamic procedure for computation of the model coefficients? It's standard practice in combustion LES, but I'm not not sure about atmospheric applications. If we're mainly planning to rely on the Deardorff model anyway, maybe it is not needed, or maybe it would also be useful for computing those model coefficients.

Deardorff
~~~~~~~~~
	  
In Deardorff (one-equation/TKE) models, the turbulent viscosity is computed as: 

.. math::

   \nu_t = C_k \Delta (k^{sfs})^{1/2}

and a prognostic (transport) equation is solved for the subfitler turbulent kinetic energy (:math:`k^{sgs}`):

.. math::

   \frac{\partial \overline{\rho} k^{sfs}}{\partial t} = - \nabla \cdot (\overline{\rho} \mathbf{\tilde{u}} \tilde{k}^{sfs}) + \nabla \cdot \left( \overline{\rho} \frac{\nu_t}{\sigma_k} \nabla k ^{sfs}  \right) + ( \overline{\rho} \widetilde{\mathbf{uu}} - \overline{\rho} \tilde{\mathbf{u}} \tilde{\mathbf{u}})\nabla \cdot \mathbf{\tilde{u}} - \overline{\rho} C_\epsilon \frac{(k^{sfs})^{3/2}}{\overline{\Delta}}

.. note:: Additional terms required for bouyancy?
   
In this modeled form of the equation, viscous transport has been neglected, turbulent transport has been modeled with an eddy viscosity model, and dissipation is modeled based on dimensional arguments. The prevgiously definted model for :math:`\tau^{sfs}_{ij}` is necessary for closing the productiont term.

.. note:: It's still not clear to me where ABL/PBL schemes fit in. It would be helpful if someone could point out where they would show up in the above equations (and/or their boundary conditions).
   
Here :math:`\rho, T, \theta`, and :math:`p` are the density, temperature, potential temperature and pressure, respectively;
these variables are all defined at cell centers.
:math:`C` is an advected quantity, i.e., a tracer, also defined at cell centers.
:math:`\mathbf{u}` and :math:`(\rho \mathbf{u})` are the velocity and momentum, respectively,
and are defined on faces. :math:`R_d` and :math:`c_p` are the gas constant and specific heat capacity for dry air respectively,
and :math:`\gamma = c_p / (c_p - R_d)` .  :math:`p_0` is a reference value for pressure. Overbars indicate filtering and tildes
indicate density-weighted (Favre) filtering (e.g., :math:`\tilde{\theta} = \overline{\rho \theta} / \overline{\rho}`).


Other Considerations
--------------------
     
The equations can be re-written in perturbational form by replacing the z-momentum equation with

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


