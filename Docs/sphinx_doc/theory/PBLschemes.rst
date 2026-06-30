
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _PBLschemes:

PBL Schemes
===========

Planetary Boundary Layer (PBL) schemes are used to model unresolved transport
in the vertical direction within the planetary boundary layer when mesh
resolutions are too coarse to resolve any of the turbulent eddies responsible
for this transport (~1 km grid resolution or larger). The PBL scheme is used to
provide closure for vertical turbulent fluxes
(i.e., :math:`\widetilde{w'\phi'} = \widetilde{w\phi} - \widetilde{w}\widetilde{\phi}`,
for any quantity :math:`\phi`). PBL schemes may be used in
conjunction with an LES model that specifies horizontal turbulent transport, in
which case the vertical component of the LES model is ignored.

Right now, ERF supports several PBL schemes: MYNN Level 2.5, MYJ, SHOC, MRF, and YSU.

The MYNN Level 2.5 model is the Mellor-Yamada-Nakanishi-Niino Level 2.5 model, largely matching the original forumulation proposed by Nakanishi and Niino in a series of papers from 2001 to 2009.

.. _MYNN25:

MYNN Level 2.5 PBL Model
------------------------
In this model, the vertical turbulent diffusivities are computed in a local
manner based on a gradient diffusion approach with coefficients computed based on
a transported turbulent kinetic energy value. The implementation and description
here largely follows `Nakanishi and Niino, Journal of Meteorological Society of Japan, 2009
<https://doi.org/10.2151/jmsj.87.895>`_, but has also been influenced by the full series
of papers that led to the development of this model and by a few documents published since then,
as listed in the :ref:`Useful References<MYNNReferences>` section below.

..
  A key difference is conservative form


The prognostic equation
for :math:`q^2 = \widetilde{u_i u_i} - \widetilde{u}_i\widetilde{u}_i` is

.. math::

    \frac{\partial \bar{\rho} q^2}{\partial t}
    + \left[ \frac{\partial \bar{\rho} \widetilde{u}_i q^2}{\partial x_i} \right]
    = \frac{\partial}{\partial z} \left(K_{q,v} \frac{\partial q^2}{\partial z} \right)
    + 2\bar{\rho} \left(-\widetilde{u'w'} \frac{\partial \widetilde{u}}{\partial z}
    - \widetilde{v'w'}\frac{\partial \widetilde{v}}{\partial z}
    + \beta g \widetilde{w'\theta'}
    - \frac{q^3}{B_1 l}
    \right)

where :math:`B_1` is a model parameter, :math:`\beta` is the thermal
expansion coefficient and `l` is a lengthscale. The vertical turbulent transport
coefficients are then computed:

.. math::

   K_{m,v} = l q S_m, K_{q,v} = 3 l q S_m, K_{\theta, v} = l q S_\theta

where :math:`S_m` and :math:`S_\theta` are stability parameters thaat
account for buoyancy effects. These
coefficients are then applied in evaluating the vertical component of turbulent
fluxes in a similar manner as is described for the
:ref:`Smagorinsky LES model<SmagorinskyModel>`. Computation of the stability parameters
and lengthscale depend on the Obukhov length and surface heat flux, which are
obtained from the :ref:`sec:MOST<sec:surface_layer>`. Further detail on these
computations can be found in the cited works. Several model coefficients are
required, with default values in ERF taken from the work of Nakanishi and Niino.

.. _MYNNReferences:

Useful References
~~~~~~~~~~~~~~~~~

The following references have informed the implementation of the MYNN PBL model in ERF:

.. _Mellor73: https://doi.org/10.1175/1520-0469(1973)030<1061:APOTPO>2.0.CO;2

.. _MY74: https://doi.org/10.1175/1520-0469(1974)031<1791:AHOTCM>2.0.CO;2

- `Mellor, Journal of the Atmospheric Sciences, 1973 <Mellor73_>`_: Introduces a PBL model based on :math:`q^2`

- `Mellor and Yamada, Journal of the Atmospheric Sciences, 1974 <MY74_>`_: Introduces PBL Model Hierarchy (Levels 1-4)

- `Mellor and Yamada, Reviews of Geophysics and Space Physics, 1982 <https://doi.org/10.1029/RG020i004p00851>`_:
  Introduces Level 2.5 Model

- `Nakanishi, Boundary-Layer Meteorology, 2001 <https://doi.org/10.1023/A:1018915827400>`_: Fits new model
  coefficients and proposes new diagnostic equation for the length scale

- `Nakanishi and Niino, Boundary-Layer Meteorology, 2004 <https://doi.org/10.1023/B:BOUN.0000020164.04146.98>`_:
  Extends the MYNN PBL modeling framework for moist conditions

- `Nakanishi and Niino, Boundary-Layer Meteorology, 2006 <https://doi.org/10.1007/s10546-005-9030-8>`_:
  Numerical stability improvements for the MYNN PBL modeling framework

- `Nakanishi and Niino, Journal of the Meteorological Society of Japan, 2009 <https://doi.org/10.2151/jmsj.87.895>`_:
  Summary of MYNN model development,
  re-evaluation of coefficients, and additional demonstration cases

- `Skamarock et al., A Description of the Advanced Research WRF Model Version 4, 2021 <http://dx.doi.org/10.5065/1dfh-6p97>`_:
  Description of the models implemented in WRF

- `Olson et al., A Description of the MYNN-EDMF Scheme and the Coupling to Other Components in WRF–ARW, 2019
  <https://doi.org/10.25923/n9wm-be49>`_:
  Description of more recent advancements upon the MYNN model

- `Juliano et al., Monthly Weather Review, 2022 <https://doi.org/10.1175/MWR-D-21-0164.1>`_:
  Description of a 3D generalization Mellor-Yamada PBL models

Discussions with Branko Kosovic (NCAR) and Joseph B. Olson (NOAA) have also played a major role in informing
the implementation of MYNN PBL models in ERF.

.. _MYNNEDMF:

MYNN-EDMF Level 2.5 PBL Model
-----------------------------

.. warning::

   Implementation is in progress with basic support.

More recent advancements that add significant complexity to the MYNN scheme have been incorporated into WRF, as described in Olson et al. 2019. These advancements are not included in ERF, but may be in the future.

.. _MYJ:

MYJ PBL Model
-------------

.. warning::

   Implementation is in progress with basic support.

The Mellor-Yamada-Janjic (MYJ) scheme is a 1.5-order turbulence closure that solves
a prognostic equation for turbulent kinetic energy (TKE). It uses a local closure approach
with no counter-gradient terms, making it particularly effective for stable and neutral boundary layers.

The turbulent fluxes are computed using gradient diffusion:

.. math::
   \overline{w'\phi'} = -K_\phi \frac{\partial \phi}{\partial z}

The vertical turbulent transport coefficients are computed from TKE and a master length scale:

.. math::
   K_{m,v} = \rho L q S_m, \quad K_{\theta,v} = \rho L q S_h, \quad K_{q,v} = \rho L q S_h

where :math:`q = \sqrt{2\cdot\text{TKE}}`, :math:`L` is the master length scale, and
:math:`S_m`, :math:`S_h` are stability functions that account for buoyancy effects and depend on
the gradient Richardson number. The master length scale :math:`L` is diagnosed based on the
PBL height, von Kármán's constant, and height above the surface within the PBL, transitioning
to a local mixing length in the free atmosphere.

The prognostic TKE equation includes production by shear and buoyancy, and dissipation:

.. math::
   \frac{\partial \text{TKE}}{\partial t} + \nabla \cdot (\mathbf{u} \text{TKE})
   = P_s + P_b - \epsilon + \nabla \cdot (K_q \nabla \text{TKE})

where :math:`P_s` is shear production, :math:`P_b` is buoyancy production, and
:math:`\epsilon` is dissipation.

Closure coefficients are taken from Janjić (2002) NCEP Office Note 437. The implementation in ERF follows Janjić (1994, 2002) and uses the Mellor-Yamada (1982) length scale formulation.

References
~~~~~~~~~~

* Janjić, Z. I. (1994): "The Step-Mountain Eta Coordinate Model: Further developments
  of the convection, viscous sublayer, and turbulence closure schemes",
  *Monthly Weather Review*, 122(5), 927-945.
* Janjić, Z. I. (2002): "Nonsingular implementation of the Mellor-Yamada Level 2.5 Scheme
  in the NCEP Meso model", NCEP Office Note No. 437.
* Mellor, G. L., & Yamada, T. (1982): "Development of a turbulence closure model for
  geophysical fluid problems", *Reviews of Geophysics*, 20(4), 851-875.

.. _SHOC:

SHOC PBL Model
--------------

.. warning::

   Implementation is in progress with basic support.

The Simplified Higher-Order Closure (SHOC) is a unified parameterization that represents
both turbulent mixing and shallow convection in a single framework. Originally developed for
the Community Atmosphere Model (CAM) and now used in E3SM, SHOC uses prognostic TKE with
diagnostic second and third-order moments and assumed probability density functions (PDFs)
to represent subgrid-scale variability.

SHOC computes vertical turbulent fluxes for momentum, heat, and moisture, along with
subgrid-scale cloud fraction and liquid water content. The assumed PDFs allow the scheme
to predict partial cloudiness and transitions between clear and cloudy conditions. The
implementation uses higher-order closure equations to diagnose eddy diffusivities and
turbulent fluxes, with special treatment for cloud-top entrainment.

References
~~~~~~~~~~

* Golaz, J.-C., et al. (2002): "A PDF-based model for boundary layer clouds. Part I:
  Method and model description", *Journal of the Atmospheric Sciences*, 59(24), 3540-3551.
* Bogenschutz, P. A., & Krueger, S. K. (2013): "A simplified PDF parameterization of
  subgrid-scale clouds and turbulence for cloud-resolving models",
  *Journal of Advances in Modeling Earth Systems*, 5(2), 195-211.
* E3SM SHOC Documentation: https://github.com/E3SM-Project/E3SM/tree/master/components/eamxx/src/physics/shoc

.. _MRFPBL:

MRF PBL Model
-------------

.. warning::

   Implementation is in progress with basic support. Need to be tuned in future for real flows.

The Medium Range Forecast (MRF) PBL model is a nonlocal boundary layer scheme that was originally
developed for the MRF model, which was used in the NCEP global forecast system. It is a nonlocal
scheme that uses a countergradient diffusion approach to model vertical turbulent transport within
the PBL. The implementation in ERF follows the original Hong and Pan (1996) formulation with several
modern enhancements described below.

**Primary Reference:** Hong, S. Y., and H.-L. Pan, 1996: Nonlocal Boundary Layer Vertical
Diffusion in a Medium-Range Forecast Model. *Monthly Weather Review*, 124, 2322-2339.
https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2

**WRF Implementation:** module_bl_mrf.F
https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F

The turbulent diffusion for prognostic variables (:math:`C= u, v, \theta, q_k`), where :math:`q_k`
includes all moisture variables is given by

.. math::
   \frac{\partial C}{\partial t}
   = \frac{\partial}{\partial z} \left[
   K_c \left( \frac{\partial C}{\partial z} - \gamma_c \right)
   \right]

Here :math:`K_c` is the turbulent diffusion coefficient, and :math:`\gamma_c` is the countergradient
correction term (nonlocal flux).

The turbulent diffusion coefficient in the mixed layer is given by:

.. math::
   K_m = \kappa w_s z \left( 1- \frac{z}{h} \right)^2

.. math::
   w_s = \frac{u_*}{\phi_m}

where :math:`\kappa` is the von Karman constant, :math:`w_s` is a representative velocity scale in
the mixed layer, and :math:`h` is the PBL height. The stability function :math:`\phi_m` is computed
to be consistent with the surface layer bottom. For unstable regime (:math:`u_*\theta_* < 0`), it is
calculated as follows:

.. math::
   \phi_m = \left(1 - 8 sf \frac{h}{L}\right)^{-1/3}

.. math::
   \phi_{t,q} = \left(1 - 16 sf \frac{h}{L}\right)^{-1/2}

and for stable regime (:math:`u_*\theta_* > 0`), it is calculated as:

.. math::
   \phi_{m,t,q} = \left(1 + 5 sf \frac{h}{L}\right)

where :math:`sf`  is a fraction of the surface layer and atmospheric boundary layer height and
:math:`L` is the Monin-Obukhov length, which is computed from the surface heat fluxes. The turbulent
coefficient for temperature and moisture is given by:

.. math::
   K_t = K_q = \frac{K_m}{Pr}

.. math::
   Pr = \left(\frac{\phi_t}{\phi_m}+ b \kappa sf\right)

where :math:`K_t` is the turbulent diffusion coefficient for temperature, :math:`K_q` is the
turbulent diffusion coefficient for moisture and :math:`Pr` is the Prandtl number.

The turbulent diffusion coefficient in the free atmosphere is computed from the YSU model as the MRF
expressions showed oscillations in the canonical stable boundary layer tests.

.. math::
   K_{m,t} = l^2 f_{m,t}(Rig)\left|\frac{\partial U}{\partial z}\right|

.. math::
   l = \frac{\kappa z \lambda}{\kappa z + \lambda}

where :math:`l` is the length scale, :math:`f_{m,t}` is a stability function for momentum and
temperature (or moisture), :math:`Rig` is the gradient Richardson number,  and :math:`U` is the
horizontal wind speed. The gradient Richardson number is computed as:

.. math::
   Rig = \frac{g}{\theta}\left[\frac{\partial \theta}{\partial z} \left(\frac{\partial U}{\partial z}\right)^{-2}\right]

A different expression is used for the stability function :math:`f_{m,t}` for stable and unstable
regimes. For stable regime we have,

.. math::
   f_t = f_m (1+2.1 Rig) = \frac{1}{\left(1 + 5 Rig\right)^2}

For the unstable regime, we have:

.. math::
   f_t = 1 - \frac{8 Rig}{1+1.286\sqrt{-Rig}}

.. math::
   f_m = 1 - \frac{8 Rig}{1+1.746\sqrt{-Rig}}

The countergradient correction term is given by:

.. math::
   \gamma_c = b \frac{ u_* \theta_*}{w_s}

where :math:`b=7.8` is a constant, :math:`u_*` is the surface frictional velocity scale, :math:`\theta_*`
is the surface potential temperature scale.

.. note::

   The countergradient correction term is now optional in ERF and can be enabled via the
   ``enable_mrf_countergradient`` flag in the input file. By default, it is disabled to maintain
   backward compatibility. Set ``erf.enable_mrf_countergradient = true`` to enable this feature.

.. _MRFEnhancements:

MRF Model Enhancements and Extensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MRF model has been enhanced with several optional features in ERF that extend beyond the original
Hong and Pan (1996) formulation, aligning more closely with WRF's advanced parameterizations and
modern understanding of boundary layer physics. All features are **disabled by default** to maintain
backward compatibility with existing simulations. This section documents all enhancements and their
physical justification.

**Important Note on Physics Fidelity:**
When using the original MRF scheme without enhancements, users should cite Hong and Pan (1996) and
acknowledge that they are using the standard formulation. When any enhancement is enabled, proper
attribution to ERF documentation and the referenced papers is required.

1. VPERT (Virtual Potential Temperature Perturbation) Correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Issue Fixed:** The original WRF MRF code (module_bl_mrf.F line 878) includes the following operation:

.. code-block:: fortran

   VPERT = HGAMT + EP1*THX(I,KL)*HGAMQ
   VPERT = MIN(VPERT, GAMCRT)  ! <-- This is incorrect physics

**Problem:** After combining sensible heat (HGAMT) and latent heat (HGAMQ) contributions, limiting
VPERT to GAMCRT suppresses the effect of latent heating. This violates the physics because:

1. HGAMT is already limited to GAMCRT (max sensible heating)
2. Adding positive HGAMQ can legitimately produce VPERT > GAMCRT
3. The latent heating from evaporation in unstable conditions is being artificially suppressed
4. This causes underprediction of PBL height in strongly moist convective conditions

**ERF Correction:** ERF can implement unbounded VPERT when ``enable_mrf_unbounded_vpert`` is enabled:

.. code-block:: cpp

   const Real VPERT = amrex::max(HGAMT + 0.61 * theta * HGAMQ, zero);

Otherwise, for strict WRF consistency, it bounds VPERT to ``GAMCRT`` (3 K):

.. code-block:: cpp

   const Real VPERT = amrex::max(amrex::min(HGAMT + 0.61 * theta * HGAMQ, GAMCRT), zero);

Enabling unbounded VPERT preserves the combined heating effect while preventing negative VPERT values (via MAX with zero).
The correction produces physically more accurate PBL heights in moist environments.

**References:**
- Original error identified during ERF development
- See ERF_ComputeDiffusivityMRF.cpp source comments for detailed physics explanation

2. HGAMQ Moisture Countergradient with Proper Limiting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Sign Convention Clarification:** The WRF convention uses:
- :math:`q_*` is negative for upward moisture flux (evaporation from surface)
- :math:`u_*` is positive (friction velocity magnitude)
- Therefore: :math:`-const_b \cdot u_* \cdot q_*` gives positive HGAMQ for unstable (evaporating) conditions

**Missing WRF Safeguard:** WRF applies MAX limiting (HGAMQ = MAX(HGAMQ, 0.0)) but only after the
MIN operation. ERF correctly applies full bounding:

.. code-block:: cpp

   HGAMQ = max(min(-const_b * u_* * q_* / w_*, GAMCRQ), zero);

This prevents unrealistic negative moisture countergradient if :math:`q_*` becomes positive (condensation),
which would indicate upside-down countergradient flux.

**Land/Water Discrimination:** HGAMQ is zeroed over water surfaces because:
- Evaporation over water is implicitly handled by ocean/water body parameterizations
- Countergradient moisture fluxes are primarily a land phenomenon (soil moisture limitation)
- This follows WRF's approach (module_bl_mrf.F line 876)

3. Cloud-Aware Stability Function Adjustments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose:** Improves representation of cloudy boundary layers by modulating stability functions
based on cloud water/ice content.

**Implementation Details:**

- Cloud detection threshold: :math:`q_c + q_i > 0.1 \text{ g/kg}` (1e-4 kg/kg)
- **In stable layers with clouds:** Reduce stability damping by 10-20% where cloud content exceeds threshold

  * Physics: Clouds reduce vertical oscillations through radiative and latent effects, making the
    layer less suppressive to turbulence
  * Improved representation of fog/stratus-topped stable layers

- **In unstable layers with clouds:** Slightly enhance instability enhancement by 5% where cloud
  content exceeds threshold

  * Physics: Latent heat release from condensation enhances buoyancy
  * Better captures cumulus-capped boundary layers

**Parameters:**
- Enabled via: ``enable_mrf_countergradient`` flag (default: false)
- Adjustment strength: 15-20% reduction in stable, 5% boost in unstable
- Can be customized via ``pbl_mrf_cloud_adjustment_factor`` parameter

**Physical Justification:**
- Clouds modify vertical buoyancy structure through radiative cooling/warming
- Latent heat release enhances convective mixing
- Cloud-top entrainment zones are qualitatively different from clear-air turbulence
- Conceptually similar to WRF's IMVDIF cloud-aware parameterization (Bretherton & Park 2009)

**References:**
- Bretherton, C. S., and S. Park, 2009: A new moist turbulence parameterization in the WRF
  Advanced Research WRF (ARW) model. In *Proceedings of the 9th Annual WRF Users' Workshop*.

4. Virtual Potential Temperature Treatment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Enhancement:** Proper handling of moisture effects on buoyancy throughout the scheme.

- Computes :math:`\theta_v = \theta(1 + 0.61 \cdot q_v)` at all levels
- Used in Richardson number diagnosis: :math:`Rig = \frac{g}{\theta_v} \frac{d\theta_v}{dz} / (shear)^2`
- Critical for accurate PBL height and convective strength

This is standard meteorological practice and implicitly assumed in Hong & Pan (1996), now made explicit.

5. Free Atmosphere Mixing via YSU Stability Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Reason for Enhancement:** Original MRF diffusivity formulation caused oscillations in canonical
stable boundary layer tests (GABLS cases).

**Solution:** Use YSU scheme (Hong et al. 2006, Appendix A) Richardson number-dependent mixing above PBL:

- More stable numerically
- Better represents free atmosphere mixing
- Standard in modern WRF configurations
- Detailed equations and references provided in source code comments

**References:**
- Hong, S. Y., Y. Noh, and J. Dudhia, 2006: A new vertical diffusion package with an explicit
  treatment of entrainment processes. *Monthly Weather Review*, 134, 2318-2341.
  https://doi.org/10.1175/MWR3250.1

Older MRF Enhancements (Deprecated/Documented for Historical Completeness)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These features are historically documented but superseded by the above enhancements:

Iterative Thermal Excess Correction
""""""""""""""""""""""""""""""""""""

**Purpose:** The thermal excess correction (:math:`\theta_T`) modifies the surface virtual potential
temperature to account for cumulative heating effects in the PBL. The iterative method refines this
estimate through multiple passes, similar to WRF's HGAMT/HGAMQ formulation, providing better estimates
than the single-pass simple method:

.. math::

   \theta_T = -b \frac{u_* \theta_*}{w_*}

**Parameters:**

- ``pbl.enable_mrf_iterative_thermal_excess`` (bool): Enable iterative refinement (default: false)
- ``pbl.mrf_thermal_excess_iterations`` (int): Number of refinement iterations (default: 3)

**Usage Example:**

.. code-block:: bash

   pbl.pbl_type = "mrf"
   pbl.enable_mrf_iterative_thermal_excess = true
   pbl.mrf_thermal_excess_iterations = 3

**Physical Justification:**

- Simple method: Valid for most stable/neutral conditions
- Iterative method: Better for strong convective conditions with multiple heating cycles
- Convergence: Typically reaches solution within 3-5 iterations

Saturated Layer Handling (IMVDIF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose:** Modulates moisture diffusivity in cloud/saturated layers using a height-dependent cloud
fraction estimate. This addresses the reduced diffusivity characteristic of stratocumulus-topped
boundary layers and fog/stratus evolution.

**Parameters:**

- ``pbl.enable_mrf_cloudy_layers`` (bool): Enable cloud-aware modulation (default: false)
- ``pbl.mrf_cloud_diffusivity_factor`` (Real): Diffusivity reduction factor in cloudy layers (default: 0.8)

  - Range: [0, 1]
  - 1.0 = no reduction (default diffusivity)
  - 0.8 = 20% reduction in cloudy regions
  - 0.0 = complete suppression

**Implementation Details:**

- Uses a simple heuristic: reduces diffusivity in lower layers (z < 2000 m)
- Gradual transition: full reduction near surface, tapering upward
- Can be enhanced with actual cloud fraction variables from microphysics

**Usage Example:**

.. code-block:: bash

   pbl.pbl_type = "mrf"
   pbl.enable_mrf_cloudy_layers = true
   pbl.mrf_cloud_diffusivity_factor = 0.8

**Physical Justification:**

- Stratocumulus layers are less turbulent than well-mixed layers
- Reduced diffusivity better captures cloud-top entrainment dynamics
- Improves representation of fog/stratus evolution

Countergradient Term Bounding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose:** Applies bounds to countergradient (non-local) turbulent flux terms, preventing unrealistic
amplification of fluxes. This follows WRF's approach with GAMCRT and GAMCRQ parameters.

**Parameters:**

- ``pbl.enable_mrf_countergradient_bounds`` (bool): Enable bounds (default: false)
- ``pbl.mrf_countergradient_max_theta`` (Real): Maximum heat countergradient (default: 3.0)

  - WRF equivalent: GAMCRT = 3

- ``pbl.mrf_countergradient_max_q`` (Real): Maximum moisture countergradient (default: 0.002)

  - WRF equivalent: GAMCRQ = 2E-3

**Implementation Details:**

The bounds are applied as:

.. code-block:: cpp

   cg_theta = min(countergradient_theta, mrf_countergradient_max_theta)
   cg_q = min(countergradient_q, mrf_countergradient_max_q)

**Usage Example:**

.. code-block:: bash

   pbl.pbl_type = "mrf"
   pbl.enable_mrf_countergradient = true
   pbl.enable_mrf_countergradient_bounds = true
   pbl.mrf_countergradient_max_theta = 3.0
   pbl.mrf_countergradient_max_q = 0.002

**Physical Justification:**

- Prevents countergradient fluxes from exceeding realistic bounds
- Critical for very strong convective conditions
- Improves model stability in extreme parameter regimes

Iterative Thermal Excess Correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose:** The thermal excess correction (:math:`\theta_T`) modifies the surface virtual potential temperature to account for cumulative heating effects in the PBL. The iterative method refines this estimate through multiple passes, similar to WRF's HGAMT/HGAMQ formulation, providing better estimates than the single-pass simple method:

.. math::

   \theta_T = -b \frac{u_* \theta_*}{w_*}

**Parameters:**

- ``pbl.enable_mrf_iterative_thermal_excess`` (bool): Enable iterative refinement (default: false)
- ``pbl.mrf_thermal_excess_iterations`` (int): Number of refinement iterations (default: 3)

**Usage Example:**

.. code-block:: bash

   pbl.pbl_type = "mrf"
   pbl.enable_mrf_iterative_thermal_excess = true
   pbl.mrf_thermal_excess_iterations = 3

**Physical Justification:**

- Simple method: Valid for most stable/neutral conditions
- Iterative method: Better for strong convective conditions with multiple heating cycles
- Convergence: Typically reaches solution within 3-5 iterations

Saturated Layer Handling (IMVDIF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose:** Modulates moisture diffusivity in cloud/saturated layers using a height-dependent cloud fraction estimate. This addresses the reduced diffusivity characteristic of stratocumulus-topped boundary layers and fog/stratus evolution.

**Parameters:**

- ``pbl.enable_mrf_cloudy_layers`` (bool): Enable cloud-aware modulation (default: false)
- ``pbl.mrf_cloud_diffusivity_factor`` (Real): Diffusivity reduction factor in cloudy layers (default: 0.8)

  - Range: [0, 1]
  - 1.0 = no reduction (default diffusivity)
  - 0.8 = 20% reduction in cloudy regions
  - 0.0 = complete suppression

**Implementation Details:**

- Uses a simple heuristic: reduces diffusivity in lower layers (z < 2000 m)
- Gradual transition: full reduction near surface, tapering upward
- Can be enhanced with actual cloud fraction variables from microphysics

**Usage Example:**

.. code-block:: bash

   pbl.pbl_type = "mrf"
   pbl.enable_mrf_cloudy_layers = true
   pbl.mrf_cloud_diffusivity_factor = 0.8

**Physical Justification:**

- Stratocumulus layers are less turbulent than well-mixed layers
- Reduced diffusivity better captures cloud-top entrainment dynamics
- Improves representation of fog/stratus evolution

Countergradient Term Bounding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose:** Applies bounds to countergradient (non-local) turbulent flux terms, preventing unrealistic amplification of fluxes. This follows WRF's approach with GAMCRT and GAMCRQ parameters.

**Parameters:**

- ``pbl.enable_mrf_countergradient_bounds`` (bool): Enable bounds (default: false)
- ``pbl.mrf_countergradient_max_theta`` (Real): Maximum heat countergradient (default: 3.0)

  - WRF equivalent: GAMCRT = 3

- ``pbl.mrf_countergradient_max_q`` (Real): Maximum moisture countergradient (default: 0.002)

  - WRF equivalent: GAMCRQ = 2E-3

**Implementation Details:**

The bounds are applied as:

.. code-block:: cpp

   cg_theta = min(countergradient_theta, mrf_countergradient_max_theta)
   cg_q = min(countergradient_q, mrf_countergradient_max_q)

**Usage Example:**

.. code-block:: bash

   pbl.pbl_type = "mrf"
   pbl.enable_mrf_countergradient = true
   pbl.enable_mrf_countergradient_bounds = true
   pbl.mrf_countergradient_max_theta = 3.0
   pbl.mrf_countergradient_max_q = 0.002

**Physical Justification:**

- Prevents countergradient fluxes from exceeding realistic bounds
- Critical for very strong convective conditions
- Improves model stability in extreme parameter regimes

High-Resolution Grid-Dependent Diffusivity Bounds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose:** Applies grid-dependent (high-resolution) bounds to diffusivity coefficients, permitting stronger and more resolved mixing in fine-resolution Large Eddy Simulations (LES) or high-resolution mesoscale runs ($\Delta z < 100\ \text{m}$) compared to conservative global forecast limits.

**Parameters:**

- ``pbl.pbl_mrf_highres_bounds`` (bool): Enable high-resolution grid-dependent bounds (default: false)

**Implementation Details:**

When enabled, the diffusivity bounds scale with local grid-spacing and density as formulated in Hong et al. (2006):

- Conservative Bounds (Default):
  :math:`K_{min} = 0.1\ \text{m}^2/\text{s}, K_{max} = 300\ \text{m}^2/\text{s}`
- High-Resolution Bounds:
  :math:`K_{min} = 0.001 \times \Delta z \times \rho, K_{max} = 1000\ \text{m}^2/\text{s}`

Combined Usage Examples
^^^^^^^^^^^^^^^^^^^^^^^

**Default Configuration (All Features Off):**

.. code-block:: bash

   # Standard MRF configuration (existing parameters)
   pbl.pbl_type = "mrf"
   pbl.pbl_mrf_Ribcr = 0.5
   pbl.pbl_mrf_const_b = 7.8
   pbl.pbl_mrf_sf = 0.1

Use this configuration for:

- Standard ABL simulations
- Backward compatibility
- Neutral/stable boundary layers

**Conservative Enhancement:**

.. code-block:: bash

   pbl.pbl_type = "mrf"
   pbl.enable_mrf_iterative_thermal_excess = true
   pbl.mrf_thermal_excess_iterations = 3

Use for: Improving strong convection without other changes

**Cloud-Aware Configuration:**

.. code-block:: bash

   pbl.pbl_type = "mrf"
   pbl.enable_mrf_cloudy_layers = true
   pbl.mrf_cloud_diffusivity_factor = 0.8

Use for: Stratocumulus-topped or fog/stratus simulations

**Full WRF-Style Configuration:**

.. code-block:: bash

   pbl.pbl_type = "mrf"
   pbl.enable_mrf_iterative_thermal_excess = true
   pbl.mrf_thermal_excess_iterations = 3
   pbl.enable_mrf_cloudy_layers = true
   pbl.mrf_cloud_diffusivity_factor = 0.8
   pbl.enable_mrf_countergradient = true
   pbl.enable_mrf_countergradient_bounds = true
   pbl.mrf_countergradient_max_theta = 3.0
   pbl.mrf_countergradient_max_q = 0.002

Use for: Maximum realism in challenging conditions

Performance Impact
^^^^^^^^^^^^^^^^^^

- **Iterative thermal excess**: ~1-5% additional computation
- **Cloud-aware moisture**: <1% overhead (simple height-based heuristic)
- **Countergradient bounds**: <1% overhead (minimal arithmetic)
- **Total impact with all enabled**: ~1-5% runtime increase

.. _MRFReferences:

Useful References
^^^^^^^^^^^^^^^^^

- Hong et al. (1996): "Nonlocal Boundary Layer Vertical Diffusion in a Medium-Range Forecast Model"
- Hong et al. (2006): "A new vertical diffusion package with an explicit treatment of entrainment processes"
- WRF Model Documentation: PBL Schemes
- Skamarock et al., A Description of the Advanced Research WRF Model Version 4, 2021 <http://dx.doi.org/10.5065/1dfh-6p97>

.. _YSUPBL:

YSU PBL Model
-------------

.. warning::

   Implementation is in progress, this option is not yet supported

The Yonsei University (YSU) PBL model is another commonly use scheme in WRF. It includes nonlocal mixing with  contergradient diffusion within
the PBL, and a local mixing treatment outside the PBL.

Turbulent diffusion for prognostic variables (:math:`C, u, v, \theta, q_k`), where :math:`q_k` includes all moisture variables and :math:`C`
any additional scalars (other terms in the equations omitted for brevity):

.. math::
   \frac{\partial C}{\partial t}
   = \frac{\partial}{\partial z} \left[
   K_c \left( \frac{\partial C}{\partial z} - \gamma_c \right)
   - \overline{\left(w'c' \right)_h} \left( \frac{z}{h} \right)^3
   \right]

.. note::

   Not applied for vertical velocity?

Where for each variable the turbulent diffusion coefficient :math:`K_c`, countergradient correction :math:`\gamma_c`,
and entrainment flux at the PBL top :math:`\overline{\left(w'c' \right)_h}` must be diagnosed for each variable.
The main controlling parameter is the PBL height :math:`h`.
Notably, a nonlocal model for turbulent diffusion is used for :math:`z \leq h`, but a local model is used for :math:`z \ge h`.

The first step is to determine the PBL height :math:`h`. This is defined as the smallest value of :math:`z` where the bulk
Richardson number equals the critical value, which is taken to be 0:

.. math::

   {\rm Rib}(z) = \frac{ g \left[ \theta_m(z) - \theta_s\right] }{\theta_{ma} U(z)^2}z

.. math::

   {\rm Rib}(h) = {\rm Rib_{cr}} = 0

where

- :math:`\theta_m` is the moist potential temperature,
- :math:`\theta_{ma}` is the value at the lowest vertical cell in a column,
- :math:`U = \sqrt{u^2 + v^2}` is the horizontal wind speed,
- :math:`\theta_s = \theta_{ma} + \theta_T` is the virtual temperature near the surface,
- :math:`\theta_T = a\frac{\overline{\left(w'\theta_m' \right)_0}}{w_{s0}}` is the excess virtual temperature near the surface,
- :math:`a` is a constant taken to be 6.8 per HND06 (matching the :math:`b` constant that appears elsewhere in the YSU model)
- :math:`\overline{\left(w'\theta_m' \right)_0}` is the surface virtual heat flux (determined from the MOST surface layer model),
- :math:`w_{s}(z) = \left(u_*^3 + 8 k w_{*b}^3z/h \right)^{1/3}` is a representative velocity scale in the mixed layer, with :math:`w_{s0} = w_s(h/2)` (note this equation matches the WRF implementation and description in H10, but differs from HND06, where :math:`\phi_m` appears in place of the constant 8),
- :math:`u_*` is the surface frictional velocity scale determined from the MOST surface layer model,
- :math:`k = 0.4` is the von Karman constant
- :math:`w_{*b} = \left[ g/\theta_{ma} \overline{\left(w'\theta_m' \right)_0} h \right]^{1/3}` for :math:`\overline{\left(w'\theta_m' \right)_0} > 0`, :math:`w_{*b} = 0` otherwise, is a convective velocity scale for moist air

In practice, an approximate value of :math:`h` is determined through a two-step process. First, :math:`\theta_T` is set to be zero
and a provisional value of :math:`h` is estimated. Then this provisional value of :math:`h` is used to compute :math:`\theta_T`,
which is in turn used to provide an improved estimate of :math:`h`, which is the value used in subsequent calculations.

.. note::

   This two step-process matches the WRF implementation, but this could be extended iteratively to reach convergence.


Countergradient corrections are computed as follows:

.. math::

   \gamma_\theta =

.. math::
   \gamma_u =

.. math::
   \gamma_v =

.. math::
   \gamma_{q_k} = \gamma_C = 0

Entrainment fluxes are computed:

.. math::
   \overline{\left(w'c' \right)_h} =

.. math::
   \overline{\left(w'c' \right)_h} =

Within the PBL (:math:`z \leq h`),

.. _YSUReferences:

Useful References
~~~~~~~~~~~~~~~~~

The following references have informed the implementation of the MRF and YSU model in ERF:

.. _HP96: https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2

- [H10] `Hong, Quarterly Journal of the Royal Meteorological Society, 2010 <https://doi.org/10.1002/qj.665>`_: Most up-to-date YSU model formulation as implemented in WRF, with revisions for stable boundary layers

- [HND06] `Hong, Noh, and Dudhia, Monthly Weather Review, 2006 <https://doi.org/10.1175/MWR3199.1>`_: Initial formulation referred to as the YSU model, adds improved entrainment formulation (relative to NCHR03) to work of TM86 and a few other modifications

- [NCHR03] `Noh, Cheon, Hong, and Raasch, Boundary-Layer Meteorology, 2003 <https://doi.org/10.1023/A:1022146015946>`_: Entrainment effects added to TM86

- [HP96] `Hong and Pan, Monthly Weather Review, 1996 <HP96_>`_: Largely an implementation and evaluation of TM86

- [TM86] `Troen and Mahrt, Boundary-Layer Meteorology, 1986 <https://doi.org/10.1007/BF00122760>`_: Initial incorporation of nonlocal counter-graident term in vertical diffusion model

- [WF18] `Wilson and Fovell, Weather and Forecasting, 2018 <https://doi.org/10.1175/WAF-D-17-0109.1>`_: Extension of YSU to handle interplay between radiation and fog, active in WRF with the ``ysu_topdown_pblmix = 1`` option

- The WRF Fortran source code for this `module <https://github.com/wrf-model/WRF/blob/a8eb846859cb39d0acfd1d3297ea9992ce66424a/phys/module_bl_ysu.F>`_ as of Dec. 2023. The ERF implementation supports the same physical models as this WRF implementation, with the exception of the ``ysu_topdown_pblmix = 1`` option from WF18, i.e. the implementation in ERF largely matches the PBL scheme described in H10.
