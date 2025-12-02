
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
obtained from the :ref:`sec:MOST`. Further detail on these
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

The Medium Range Forecast (MRF) PBL model is a nonlocal PBL scheme that was originally developed for the MRF model,
which was used in the NCEP global forecast system. It is a nonlocal scheme that uses a countergradient diffusion approach
to model vertical turbulent transport within the PBL.

The turbulent diffusion for prognostic variables (:math:`C= u, v, \theta, q_k`), where :math:`q_k` includes all moisture
variables is given by

.. math::
   \frac{\partial C}{\partial t}
   = \frac{\partial}{\partial z} \left[
   K_c \left( \frac{\partial C}{\partial z} - \gamma_c \right)
   \right]

Here :math:`K_c` is the turbulent diffusion coefficient, and :math:`\gamma_c` is the countergradient correction term.

The turbulent diffusion coefficient in the mixed layer is given by:

.. math::
   K_m = \kappa w_s z \left( 1- \frac{z}{h} \right)^2

.. math::
   w_s = \frac{u_*}{\phi_m}

where :math:`\kappa` is the von Karman constant, :math:`w_s` is a representative velocity scale in the mixed layer,
and :math:`h` is the PBL height. The stability function :math:`\phi_m` is computed to be consistent with the surface layer
bottom. For unstable regime (:math:`u_*\theta_* < 0`), it is calculated as follows:

.. math::
   \phi_m = \left(1 - 8 sf \frac{h}{L}\right)^{-1/3}

.. math::
   \phi_{t,q} = \left(1 - 16 sf \frac{h}{L}\right)^{-1/2}

and for stable regime (:math:`u_*\theta_* > 0`), it is calculated as:

.. math::
   \phi_{m,t,q} = \left(1 + 5 sf \frac{h}{L}\right)

where :math:`sf`  is a fraction of the surface layer and  atmospheric boundary layer height and  :math:`L`
is the Monin-Obukhov length,  which is computed from the surface heat fluxes. The turbulent coefficient for
temperature and moisture is given by:

.. math::
   K_t = K_q = \frac{K_m}{Pr}

.. math::
   Pr = \left(\frac{\phi_t}{\phi_m}+ b \kappa sf\right)

where :math:`K_t` is the turbulent diffusion coefficient for temperature, :math:`K_q` is the turbulent diffusion coefficient for moisture
and :math:`Pr` is the Prandtl number.

The turbulent diffusion coefficient in the free atmosphere is computed from the YSU model as the MRF
expressions showed oscillations in the canonical stable boundary layer tests.

.. math::
   K_{m,t} = l^2 f_{m,t}(Rig)\left|\frac{\partial U}{\partial z}\right|

.. math::
   l = \frac{\kappa z \lambda}{\kappa z + \lambda}

where :math:`l` is the length scale, :math:`f_{m,t}` is a stability function for momentum and temperature (or moisture),
:math:`Rig` is the gradient Richardson number,  and :math:`U` is the horizontal wind speed. The gradient Richardson
number is computed as:

.. math::
   Rig = \frac{g}{\theta_v}\left[\frac{\partial \theta_v}{\partial z} \left(\frac{\partial z}{\partial U}\right)^2\right]

A different expression is used for the stability function :math:`f_{m,t}` for stable and unstable regimes. For stable regime we have,

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

where :math:`b=7.8` is a constant, :math:`u_*` is the surface frictional velocity scale, :math:`\theta_*` is the
surface potential temperature scale.


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
