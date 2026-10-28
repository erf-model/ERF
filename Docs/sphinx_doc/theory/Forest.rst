.. _sec:Forest:

Forest Model
--------------

The forest model provides an option to include the drag from forested regions to be included in the momentum equation. The
drag force is calculated as follows:

.. math::

   F_i= - C_d L(x,y,z) U_i | U_i |


:math:`F_i` is an acceleration; the momentum equations, which are written for
:math:`\rho U_i`, receive :math:`\rho F_i` with the density averaged to the
face of each velocity component.

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
  ``tree_type x_center y_center height diameter Cd LAI LAImax``. Where patches
  overlap, their drag :math:`C_d L` and leaf area densities :math:`L` add. A
  cell receives drag from a patch when its centre is strictly below the tree
  height :math:`h`, as in the gridded mode; the type-2 profile tends to zero at
  :math:`z = h`.
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

Canopy turbulence in the k-eqn RANS model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The canopy drag removes mean kinetic energy from the flow. Part of it becomes
turbulence in the wakes of leaves and branches, at scales far below the mesh,
and dissipates quickly (the spectral short cut). With the one-equation
:math:`k` model (``erf.rans_type = kEqn``), ``erf.forest_canopy_tke = true``
adds both effects where the forest drag field :math:`f = C_d L` is positive:

.. math::

   S_k = \rho f \left( \beta_p |\boldsymbol{U}|^3 - \beta_d |\boldsymbol{U}| k \right),

where :math:`|\boldsymbol{U}|` is the cell-centred wind speed, :math:`\beta_p`
(``erf.forest_canopy_beta_p``, default 1) is the fraction of the drag work that
becomes TKE and :math:`\beta_d` (``erf.forest_canopy_beta_d``, default 4) sets
the short-circuit dissipation (Green 1992; Liu et al. 1996; Sanz 2003).

**Time integration.** The sink relaxes :math:`k` at the rate
:math:`c = \beta_d f |\boldsymbol{U}|`, which in a dense canopy is fast compared
with the time step. Each RK stage updates :math:`\rho k` from its value at the
start of the step, :math:`(\rho k)^n`, over :math:`\Delta t_s` (:math:`\Delta t/3`,
:math:`\Delta t/2` and :math:`\Delta t` for RK3; :math:`\Delta t` for both stages of
the anelastic RK2). With :math:`P = \beta_p f |\boldsymbol{U}|^3` and :math:`c`
taken from the stage, production and sink are integrated together and exactly
over the stage,

.. math::

   S_k = \rho P \, \varphi_1(x) - (\rho k)^n \, \frac{1 - e^{-x}}{\Delta t_s},
   \qquad x = c \, \Delta t_s, \quad \varphi_1(x) = \frac{1 - e^{-x}}{x},

so that, taken alone, the stage update is
:math:`(\rho k)^n e^{-x} + \rho (P/c)(1 - e^{-x})`: the canopy terms balance at
:math:`k = P/c` at any time step, and the canopy sink never removes more than
:math:`(\rho k)^n`. For :math:`x \ll 1` the source tends to
:math:`\rho (P - c k)`.

**Scope of the exact stage solution.** The closure's own dissipation
:math:`\rho C_{\mu 0}^3 k^{3/2}/L` is a separate term. The required implicit
form (below) divides the stage update by
:math:`1 + \Delta t_s C_{\mu 0}^3 \sqrt{k} / L`, and with
``erf.terrain_type = MovingFittedMesh`` the update is also scaled by the old and
new cell volumes. In both cases the stage value is not the expression above, and
the steady :math:`k` is the balance of all the terms rather than :math:`P/c`.
What remains true is that the canopy part of the update never makes :math:`k`
negative, and the implicit dissipation keeps it non-negative. The source is
not added in fully immersed cells (``terrain_blank`` = 1) of immersed-forcing
terrain or buildings.

Inside a dense canopy the eddies are limited by the foliage rather than by the
distance to the ground. The turbulent length scale :math:`L` of the closure
(:ref:`RANS`) is therefore capped by the canopy drag length
:math:`L_c = 1/f`,

.. math::

   L \leftarrow \min \left( L, \ \alpha L_c \right)
   \qquad \text{where } f > 0,

with :math:`\alpha` = ``erf.forest_canopy_length_alpha`` (default 0.04; zero
disables the cap and keeps the TKE source). Harman and Finnigan (2007) give
:math:`\alpha = 2 \beta^3` with :math:`\beta = u_*/U_h`, the ratio of the
friction velocity to the wind speed at the canopy top; 0.04 corresponds to
:math:`\beta \approx 0.27` and 0.054 to the often quoted :math:`\beta = 0.3`.
The cap is applied before the dissipation :math:`\rho C_{\mu 0}^3 k^{3/2}/L`,
the Richardson number and the stability functions are computed from
:math:`L`, so all of them use the capped value.

**Required: implicit TKE dissipation.** The cap makes the dissipation rate
:math:`C_{\mu 0}^3 \sqrt{k} / L = C_{\mu 0}^3 \sqrt{k} f / \alpha` large: about
1.2 1/s for :math:`f = 0.2` 1/m and :math:`k = 2` m\ :sup:`2`/s\ :sup:`2`, and about
6 1/s for :math:`f = 0.6` 1/m and :math:`k = 5` m\ :sup:`2`/s\ :sup:`2`. With the
explicit dissipation, :math:`k` oscillates and can go negative once the time step
exceeds about the inverse of that rate (0.8 s and 0.17 s in these examples).
With a nonzero :math:`\alpha`, the run therefore stops unless
``erf.implicit_tke_dissipation = true`` on every level; the implicit form has no
such limit.

The capped length has no lower bound. ``erf.rans_lscale_min`` bounds only the
PBL-height cap of the geometric length, and the closure's own length is not
bounded below either (it is :math:`\kappa (d + z_0)` next to the wall and tends to
:math:`C_b \sqrt{k}/N` in strong stratification). Typical forests have
:math:`f` = 0.02 to 0.2 1/m and so :math:`\alpha/f` = 0.2 to 2 m, which a floor of
1 m (the default of ``erf.rans_lscale_min``) would largely remove. In a very
dense canopy (:math:`f` = 1 to 2 1/m, :math:`\alpha/f` = 2 to 4 cm) the large
dissipation rate is handled by the implicit form, and the short length moves the
Richardson number towards neutral.

With ``erf.forest_substep = true`` the momentum drag moves to the acoustic
substeps, while the TKE source stays in the slow RK stages with the other
:math:`k` sources.

The input names, defaults and formulas follow the ``ForestDrag.canopy_*``
inputs of the ``KLAxell`` model in kynema-sgf (branch ``f/winsent`` of
``rthedin/amr-wind``), and in neutral air the two codes run the same canopy
model. Two details differ:

- **Stratified closure.** kynema-sgf rescales its eddy viscosity and production
  terms by the ratio of the capped to the uncapped length after they are
  computed, with the stability functions taken from the length of the previous
  step; ERF caps :math:`L` before any of them is evaluated, as in Axell and
  Liungman (2001). In neutral air the eddy viscosities agree. For a geometric
  length of 10 m, :math:`k = 0.5` m\ :sup:`2`/s\ :sup:`2` and capped lengths of
  0.2 to 0.8 m, the ERF eddy viscosity is 1.2 to 1.3 times the kynema-sgf value
  in stable air (:math:`N^2` = 1e-3 to 1e-2 1/s\ :sup:`2`; 1.7 to 2.1 times for the
  heat diffusivity) and 0.4 to 0.5 times in unstable air (0.2 to 0.3 times for the
  heat diffusivity).
- **Time integration of the TKE source.** kynema-sgf adds the production
  explicitly and integrates the sink exactly from the current :math:`k` over the
  full step; with the source evaluated once per step this settles at
  :math:`(P/c) \, x/(1 - e^{-x}) \approx (1 + c \Delta t / 2) \, P/c`, while ERF
  settles at :math:`P/c` (above). The difference is small when
  :math:`c \Delta t \ll 1`.

The canopy terms are off by default. They require a forest and
``erf.rans_type = kEqn`` on every level, and the run stops otherwise. The
coefficients below are ignored, with a warning, while
``erf.forest_canopy_tke = false``, so a deck can switch the terms on and off
with that one input.

References:

- Green, S. R. (1992). Modelling turbulent air flow in a stand of widely-spaced
  trees. PHOENICS Journal of Computational Fluid Dynamics and Its Applications,
  5, 294-312.
- Harman, I. N., and Finnigan, J. J. (2007). A simple unified theory for flow in
  the canopy and roughness sublayer. Boundary-Layer Meteorology, 123, 339-363.
- Liu, J., Chen, J. M., Black, T. A., and Novak, M. D. (1996). E-epsilon
  modelling of turbulent air flow downwind of a model forest edge.
  Boundary-Layer Meteorology, 77, 21-44.
- Sanz, C. (2003). A note on k-epsilon modelling of vegetation canopy air-flows.
  Boundary-Layer Meteorology, 108, 191-197.

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
