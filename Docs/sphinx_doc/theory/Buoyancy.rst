.. role:: cpp(code)
   :language: c++

.. role:: f(code)
   :language: fortran

.. _Buoyancy:

Buoyancy
========

ERF adds a vertical perturbational gravitational **force density** to the
face-centered dry-momentum equation. Its SI units are N m\ :sup:`-3`; it
is not an acceleration until the appropriate inertia is accounted for.
Gravity is :math:`\boldsymbol{g}=(0,0,g_z)`, where :math:`g_z=-g` and
:math:`g>0`. Positive :math:`B_z` acts upward.

The compressible runtime parameter ``erf.buoyancy_type`` accepts 1, 2, 3,
or 4. In the anelastic formulation ERF chooses a separate dry or moist
kernel, regardless of the stored numeric selector; the stored value 3 is
only bookkeeping and there is no selectable type 5. See
:ref:`sec:Inputs` and :ref:`GoverningEquations`.

Density, reference state, and water loading
-------------------------------------------

ERF predicts the **dry-air** density :math:`\rho_d` and dry potential
temperature :math:`\theta_d`. Water mixing ratios are masses per mass of
dry air. Let :math:`q_v` be water vapor, :math:`q_t` the sum of vapor and
all condensed and precipitating *water mass* mixing ratios, and
:math:`q_\ell=q_t-q_v` the total condensate and precipitate. Number
concentrations do not contribute to :math:`q_t`. The total moist density
and its hydrostatic reference are

.. math::

   \rho_t = \rho_d(1+q_t),\qquad
   \rho_{t0} = \rho_0(1+q_{v0}),\qquad
   \frac{dp_0}{dz}=\rho_{t0}g_z,

where :math:`\rho_0` is the **dry** base-state density, :math:`q_{v0}`
is base-state vapor mixing ratio, and the base state has no condensate.
After subtracting hydrostatic balance, the vertical pressure and gravity
contribution is :math:`-\partial_z(p-p_0)+B_z`. In the implemented
dry-momentum source the combined pressure and buoyancy contribution is
divided by :math:`1+q_{t,\mathrm{face}}`, where :math:`q_{t,\mathrm{face}}`
is the arithmetic mean of the two adjacent cell values.

Type 1: density perturbation
----------------------------

The density-perturbation force density is

.. math::

   B_z^{(1)} = g_z\left[\rho_d(1+q_t)-\rho_0(1+q_{v0})\right].

ERF interpolates the bracketed density perturbation from the two adjacent
cells to an interior vertical-velocity face. This formulation retains
the full compressible density response, including pressure-related
density changes through the state and equation of state. It vanishes in
the neutral moist or dry reference state. Type 1 is the default
compressible choice and is required by certain moisture-model
combinations.

Types 2 and 3: temperature-perturbation approximation
-----------------------------------------------------

**Types 2 and 3 select the same implemented kernel in dry and moist
compressible configurations.** Define

.. math::

   \epsilon_v = \frac{R_v}{R_d}-1,\qquad
   \beta_T = \frac{T-T_0}{T_0}
             +\epsilon_v(q_v-q_{v0})-q_\ell.

For dry flow the water terms vanish. The interior-face approximation is

.. math::

   B_z^{(2)}=B_z^{(3)}
      =-\overline{\rho_0}\,g_z\,\overline{\beta_T},

where each overbar is the arithmetic mean of the two adjacent cells and
:math:`T_0` is base-state temperature. This is a low-pressure-perturbation
and dilute-moisture approximation to the full density force, not an exact
replacement for type 1 in arbitrary compressible states. The moisture
term uses the vapor perturbation :math:`q_v-q_{v0}` and all
condensed/precipitating water mass :math:`q_\ell`.

Type 4: potential-temperature-perturbation approximation
--------------------------------------------------------

Define the dry potential-temperature perturbation
:math:`\theta_d'=\theta_d-\theta_0` and

.. math::

   \beta_\theta = \frac{\theta_d'}{\theta_0}
           +\epsilon_v(q_v-q_{v0})-q_\ell.

For an interior non-EB face, type 4 computes

.. math::

   B_z^{(4)}=-\overline{\rho_0}\,g_z\,\overline{\beta_\theta}.

In dry flow the moisture terms vanish. This is a distinct approximation
from the temperature formulation: it omits the explicit density
contribution of pressure perturbations and is most appropriate when
those contributions are negligible.

Relation to the compressible equation of state
----------------------------------------------

Write :math:`\alpha=R_v/R_d`, and distinguish the constant EOS reference
pressure :math:`P_{00}=10^5\,\mathrm{Pa}` from the height-dependent
hydrostatic pressure :math:`p_0(z)`. At the default internally consistent
thermodynamic constants, :math:`\Gamma=1.4` and
:math:`\kappa=R_d/C_{p,d}=(\Gamma-1)/\Gamma`, with
:math:`C_{p,d}=1004.5\,\mathrm{J\,kg^{-1}\,K^{-1}}`. The EOS is

.. math::

   p=P_{00}\left[\frac{R_d\rho_d\theta_d(1+\alpha q_v)}{P_{00}}\right]^\Gamma.

For base and perturbed states satisfying this EOS, the exact density
ratio is

.. math::

   \frac{\rho_t}{\rho_{t0}}=
   \left(\frac{p}{p_0}\right)^{1/\Gamma}
   \frac{\theta_0}{\theta_d}
   \frac{1+q_t}{1+q_{v0}}
   \frac{1+\alpha q_{v0}}{1+\alpha q_v}.

Linearizing about the base state, with :math:`\delta q_v=q_v-q_{v0}`,
gives

.. math::

   \frac{\delta\rho_t}{\rho_{t0}} =
   \frac{1}{\Gamma}\frac{p'}{p_0}
   -\frac{\theta_d'}{\theta_0}
   +\left[\frac{1}{1+q_{v0}}-\frac{\alpha}{1+\alpha q_{v0}}\right]\delta q_v
   +\frac{q_\ell}{1+q_{v0}}+O(\delta^2).

Neglecting the pressure-density term and taking the dilute-moisture limit
gives the type-4 expression above. Finite base-state humidity introduces
additional coefficient differences; the active formula is **not** the
exact finite-humidity linearization. Types 2 and 3 use :math:`T'/T_0`
instead of :math:`\theta_d'/\theta_0`; since
:math:`T'/T_0\simeq\theta_d'/\theta_0+\kappa p'/p_0`, they can differ in
sign from type 1 for a pressure-only perturbation. These are differences
in approximate gravitational source terms, not by themselves
comparisons of the complete momentum equation.

Anelastic buoyancy
------------------

Anelastic dry flow uses :math:`\rho_d=\rho_0` in the base-density
approximation and computes

.. math::

   B_{z,\mathrm{dry}}^A =
   -\overline{\rho_0}g_z
   \frac{\overline{\theta_d}-\overline{\theta_0}}
        {\overline{\theta_0}}.

Anelastic moist flow uses the active potential-temperature-perturbation
kernel

.. math::

   B_{z,\mathrm{moist}}^A =
   -\overline{\rho_0}g_z\,\overline{\beta_\theta}.

Here the bar on :math:`\beta_\theta` averages cell-wise fractional
perturbations. The dry and moist paths can differ by a small
face-interpolation truncation term on a stratified background even in
the zero-moisture limit. The anelastic pressure-gradient field is
maintained by the anelastic projection; it is not recomputed from the
compressible EOS perturbational pressure.

Applicability and references
----------------------------

Embedded-boundary buoyancy has specialized implementation constraints:
the established branches are dry anelastic or dry compressible type 1.
Some moisture models also require compressible type 1. See
``Source/SourceTerms/ERF_MakeBuoyancy.cpp`` and :ref:`sec:Inputs` for the
current restrictions. The stored selector in anelastic mode does not
request a compressible buoyancy kernel.

The temperature and potential-temperature approximations are related to
formulations discussed by `Khairoutdinov and Randall (2003)
<khairoutdinov2003cloud>`_, *Journal of the Atmospheric Sciences*, **60**,
607--625, and `Bryan and Fritsch (2002) <bryan2002benchmark>`_. Their
applicability is governed by the assumptions stated above, not by an
assertion that all four input selections are physically interchangeable.

.. _khairoutdinov2003cloud: https://journals.ametsoc.org/view/journals/atsc/60/4/1520-0469_2003_060_0607_crmota_2.0.co_2.xml
.. _bryan2002benchmark: https://journals.ametsoc.org/view/journals/mwre/130/12/1520-0493_2002_130_2917_absfmn_2.0.co_2.xml
