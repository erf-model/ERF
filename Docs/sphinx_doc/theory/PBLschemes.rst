
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
for this transport (1 km grid resolution or larger). The PBL scheme is used to
provide closure for vertical turbulent fluxes
(i.e., :math:`\widetilde{w'\phi'} = \widetilde{w\phi} - \widetilde{w}\widetilde{\phi}`,
for any quantity :math:`\phi`). PBL schemes may be used in
conjunction with an LES model that specifies horizontal turbulent transport, in
which case the vertical component of the LES model is ignored.

Right now, the only PBL scheme supported in ERF is the
Mellor-Yamada-Nakanishi-Niino Level 2.5 model.

MYNN Level 2.5 PBL Model
------------------------
In this model, the vertical turbulent diffusivities are computed in a local
manner based on a transported turbulent kinetic energy value. The model was
proposed by `Nakanishi and Niino <https://doi.org/10.1023/B:BOUN.0000020164.04146.98>`_, building on the work of `Mellor and Yamada <https://doi.org/10.1029/RG020i004p00851>`_.

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
obtained from the :ref:`MOST module<MostBoundary>`. Further detail on these
computations can be found in the cited works. Several model coefficients are
required, with default values in ERF taken from the work of Nakanishi and Niino.
