
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

Right now, the only PBL scheme supported in ERF is the Mellor-Yamada-Nakanishi-Niino
Level 2.5 model, largely matching the original forumulation proposed by Nakanishi and
Niino in a series of papers from 2001 to 2009.

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
account for bouyancy effects. These
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

More recent advancements that add significant complexity to the MYNN scheme have been incorporated into WRF, as described in Olson et al. 2019. These advancements are not included in ERF, but may be in the future.

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
- :math:`a` is a constant taken to be 6.8 per HND06 (matching the :math:`b` constant that appears elsehwere in the YSU model)
- :math:`\overline{\left(w'\theta_m' \right)_0}` is the surface virtual heat flux (determined from the MOST surface layer model),
- :math:`w_{s}(z) = \left(u_*^3 + 8 k w_{*b}^3z/h \right)^{1/3}` is a representative velocity scale in the mixed layer, with :math:`w_{s0} = w_s(h/2)` (note this equation matches the WRF implementation and description in H10, but differs from HND06, where :math:`\phi_m` appears in place of the constant 8),
- :math:`u_*` is the surface frictional velocity scale determined from the MOST surface layer model,
- :math:`k = 0.4` is the von Karman constant
- :math:`w_{*b} = \left[ g/\theta_{ma} \overline{\left(w'\theta_m' \right)_0} h \right]^{1/3}` for :math:`\overline{\left(w'\theta_m' \right)_0} > 0`, :math:`w_{*b} = 0` otherwise, is a convective velcoity scale for moist air

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

The following references have informed the implementation of the YSU model in ERF:

.. _HP96: https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2

- [H10] `Hong, Quarterly Journal of the Royal Meteorological Society, 2010 <https://doi.org/10.1002/qj.665>`_: Most up-to-date YSU model formulation as implemented in WRF, with revisions for stable boundary layers

- [HND06] `Hong, Noh, and Dudhia, Monthly Weather Review, 2006 <https://doi.org/10.1175/MWR3199.1>`_: Initial formulation referred to as the YSU model, adds improved entrainment formulation (relative to NCHR03) to work of TM86 and a few other modifications

- [NCHR03] `Noh, Cheon, Hong, and Raasch, Boundary-Layer Meteorology, 2003 <https://doi.org/10.1023/A:1022146015946>`_: Entrainment effects added to TM86

- [HP96] `Hong and Pan, Monthly Weather Review, 1996 <HP96_>`_: Largely an implementation and evluation of TM86

- [TM86] `Troen and Mahrt, Boundary-Layer Meteorology, 1986 <https://doi.org/10.1007/BF00122760>`_: Initial incorporation of nonlocal counter-graident term in vertical diffusion model

- [WF18] `Wilson and Fovell, Weather and Forecasting, 2018 <https://doi.org/10.1175/WAF-D-17-0109.1>`_: Extension of YSU to handle interplay between radiation and fog, active in WRF with the ``ysu_topdown_pblmix = 1`` option

- The WRF Fortran source code for this `module <https://github.com/wrf-model/WRF/blob/a8eb846859cb39d0acfd1d3297ea9992ce66424a/phys/module_bl_ysu.F>`_ as of Dec. 2023. The ERF implementation supports the same physical models as this WRF implementation, with the exception of the ``ysu_topdown_pblmix = 1`` option from WF18, i.e. the implementation in ERF largely matches the PBL scheme described in H10.
