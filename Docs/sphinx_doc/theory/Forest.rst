.. _sec:Forest:

Forest Model
--------------

The forest model provides an option to include the drag from forested regions to be included in the momentum equation. The
drag force is calculated as follows:

.. math::

   F_i= - C_d L(x,y,z) U_i | U_i |


Here :math:`C_d` is the coefficient of drag for the forested region and :math:`L(x,y,z)` is the leaf area density (LAD) for the
forested region. A three-dimensional model for the LAD is usually unavailable and is also cumbersome to use if there are thousands
of trees. Two different models are available as an alternative:

.. math::
   L=\frac{LAI}{h}

.. math::
   L(z)=L_m \left(\frac{h - z_m}{h - z}\right)^n  exp\left[n \left(1 -\frac{h - z_m}{h - z}\right )\right]

Here :math:`LAI` is the leaf area index and is available from measurements, :math:`h` is the height of the tree, :math:`z_m` is the location
of the maximum LAD, :math:`L_m` is the maximum value of LAD at :math:`z_m` and :math:`n` is a model constant with values  6 (below :math:`z_m`) and 0.5
(above :math:`z_m`), respectively. :math:`L_m` is computed by integrating the following equation (see `Lalic and Mihailovic (2004)
<https://doi.org/10.1175/1520-0450(2004)043%3C0641:AERDLD%3E2.0.CO;2>`__):

.. math::
   LAI = \int_{0}^{h} L(z) dz

The simplified model with uniform LAD is recommended for forested regions with no knowledge of the individual trees. LAI values can be used from
climate model look-up tables for different regions around the world if no local remote sensing data is available.

Canopy input modes
~~~~~~~~~~~~~~~~~~

ERF supports two mutually exclusive canopy descriptions:

* ``erf.forest_file`` describes discrete circular forest patches. Each row is
  ``tree_type x_center y_center height diameter Cd LAI LAImax``.
* Gridded NetCDF mode uses ``erf.forest_lai_file`` and
  ``erf.forest_height_file`` together with either ``erf.forest_cd_file`` or a
  constant ``erf.forest_cd``. The corresponding NetCDF variables are ``LAI``,
  ``height``, and ``cd``. The horizontal fields are bilinearly interpolated to
  ERF cell centers. Their horizontal dimensions, spacing, and origins must
  match. Coordinate arrays must be projected Cartesian ``x`` and ``y`` arrays
  with at least two finite, strictly increasing, uniformly spaced points in
  each direction; raw ``lon``/``lat`` arrays are rejected because no projection
  is applied. Inconsistent fields are rejected before interpolation. ERF emits
  one aggregated warning per level when the target domain extends beyond the
  source-grid footprint.

``erf.forest_tree_type = 1`` selects the uniform LAD profile, while type 2
selects the Lalic--Mihailovic profile. In gridded mode,
``erf.forest_laimax`` gives :math:`z_m/h` and must satisfy
:math:`0 \leq z_m/h < 1`. Gridded input requires ERF to be
built with NetCDF support. See :ref:`sec:Inputs` for the complete parameter
contract.

Fixed-leaf-temperature exchange
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optional ``erf.forest_biophysics_heat`` path is deliberately scoped. Enable
it together with ``erf.forest_biophysics`` and prescribe leaf potential
temperature with ``erf.forest_leaf_theta_fixed``. The resulting leaf--air
temperature difference drives canopy sensible heat and, when moisture is
active, transpiration source terms. This path does **not** solve a leaf energy
balance or include radiation, photosynthesis, or prognostic stomatal behavior.

At each canopy cell, actual leaf and air temperatures are

.. math::

   T_l = \theta_l \Pi, \qquad T_a = \theta_a \Pi,

where :math:`\Pi` is the base-state Exner function. The leaf boundary resistance
uses local cell-centered wind speed,

.. math::

   r_b = 100 \sqrt{\frac{0.05}{\max(|\boldsymbol{u}|,0.01)}}.

The per-leaf-area sensible heat flux and its conserved potential-temperature
source are

.. math::

   H_l = \frac{2\rho c_p(T_l-T_a)}{r_b}, \qquad
   \frac{\partial(\rho\theta)}{\partial t}\bigg|_{canopy}
   = \frac{L H_l}{c_p\Pi}.

For moist configurations, the fixed minimum stomatal conductance
:math:`g_s=0.01\ \mathrm{mol\,m^{-2}\,s^{-1}}` defines
:math:`r_s=p/(g_s R T_a)`. The latent heat flux and vapor source are

.. math::

   LE_l = \frac{\rho L_v\max(q_{sat}(T_l,p)-q_v,0)}{1.075r_b+r_s}, \qquad
   \frac{\partial(\rho q_v)}{\partial t}\bigg|_{canopy}
   = \frac{L LE_l}{L_v}.

Here :math:`L` is LAD. The one-sided vapor deficit means this path represents
transpiration only; it does not represent dew deposition on leaves.
