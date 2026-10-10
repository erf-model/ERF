.. _Buoyancy:

Buoyancy
========

Buoyancy describes how a parcel's density differs from that of its
surroundings and how gravity acts on that difference. In ERF, the buoyancy
term is the **perturbational gravitational force per unit volume** in the
vertical momentum equation. A positive value acts upward; a negative value
acts downward. It is not, by itself, the parcel's vertical acceleration:
pressure gradients and other momentum sources also contribute.

ERF offers several definitions of this force in **compressible** simulations.
Some use the full density departure from a hydrostatic reference state;
others approximate that departure using temperature or potential
temperature. The **anelastic** equations use separate buoyancy expressions.

Choosing a formulation
----------------------

In a compressible simulation, ``erf.buoyancy_type`` controls the expression
used for the gravitational source. Gravity must also be enabled with
``erf.use_gravity = true``; selecting a buoyancy type does not turn on
gravity. With gravity disabled, the buoyancy force is zero.

.. list-table:: Compressible buoyancy options
   :header-rows: 1
   :widths: 12 35 53

   * - Value
     - Formulation
     - Scientific interpretation
   * - ``1`` (default)
     - Total-density perturbation
     - Retains the full density departure, including water loading and
       pressure-related density changes through the compressible equation of
       state. This is ERF's default compressible option.
   * - ``2`` or ``3``
     - Temperature perturbation
     - The two values invoke **the same implementation**, in both dry and
       moist compressible runs. This is an approximate buoyancy expression,
       not the full compressible density force.
   * - ``4``
     - Potential-temperature perturbation
     - Uses dry potential-temperature anomalies, vapor anomalies, and
       condensate loading. It neglects an explicit pressure-density
       contribution and is therefore an approximation.

The choice is not merely a change in notation: these expressions can give
different gravitational source terms for the same compressible state.
Unless an approximation is part of the intended model formulation, the
default type 1 retains the more complete density response.

For ``erf.anelastic = 1``, ERF selects its dry or moist **anelastic** kernel
irrespective of ``erf.buoyancy_type``. Internally, the option is stored as
``3`` at anelastic levels, but this does **not** select the compressible
type-3 formula. The input accepts values 1--4 only; there is no type 5.
See :ref:`sec:Inputs` for related configuration options.

Reference state, moisture, and sign convention
----------------------------------------------

ERF uses **dry-air** density :math:`\rho_d` and **dry** potential temperature
:math:`\theta_d`. Water mixing ratios are measured per unit mass of *dry*
air, not per unit mass of moist air. Write

.. math::

   q_t = q_v + q_{\mathrm{cond}},\qquad
   \rho_t = \rho_d(1+q_t),

where :math:`q_v` is water-vapor mixing ratio, :math:`q_t` is the total
water-mass mixing ratio supplied to the buoyancy calculation, and
:math:`q_{\mathrm{cond}}=q_t-q_v` includes represented cloud and precipitation mass in both liquid and ice
phases. This last quantity is **not just
liquid cloud water**. ERF forms :math:`q_t` from the moisture model's
conserved water-mass components; number concentrations and unrelated
non-water species are not included in that sum.

Let :math:`\rho_0(z)`, :math:`\theta_0(z)`, :math:`p_0(z)`, and
:math:`q_{v0}(z)` denote the reference dry density, potential temperature,
pressure, and vapor mixing ratio. There is no reference condensate, so the
reference *total moist* density is

.. math::

   \rho_{t0}=\rho_0(1+q_{v0}).

For a hydrostatic reference state,

.. math::

   \frac{dp_0}{dz}=-\rho_{t0}g,

where :math:`g>0` is the magnitude of gravity. We take height :math:`z`
positive upward and write :math:`g_z=-g`. Thus a negative density
perturbation gives a **positive (upward)** gravitational force.

ERF evaluates the vertical buoyancy term at vertical-velocity faces.
Throughout this page, an overbar denotes the arithmetic average of values
from the two adjacent cell centers; the formulas below refer to interior,
non-embedded-boundary faces unless stated otherwise.

The force and the momentum tendency
-----------------------------------

The buoyancy term :math:`B_z` has units of
:math:`\mathrm{N\,m^{-3}}=\mathrm{kg\,m^{-2}\,s^{-2}}`, rather than
:math:`\mathrm{m\,s^{-2}}`. On a flat mesh, omitting other forces and any
additional prescribed pressure forcing, the contribution of pressure and
buoyancy to the dry vertical-momentum tendency is

.. math::

   \left.\frac{\partial(\rho_* w)}{\partial t}\right|_{p,B}
   =\frac{-G_{p,z}+B_z}{1+q_{t,f}},\qquad
   q_{t,f}=\frac{q_t(k-1)+q_t(k)}{2},

where :math:`\rho_*` denotes :math:`\rho_d` in compressible mode and the
fixed reference dry density :math:`\rho_0` in anelastic mode. The quantity
:math:`G_{p,z}` is ERF's vertical pressure-gradient term. For compressible
flow it is calculated from the EOS pressure perturbation
:math:`p'=p-p_0`; for anelastic flow it is maintained through the
pressure projection. Anelastic projection also applies a separate
momentum correction, which is not represented by this isolated slow-RHS
expression. Terrain metrics and other forcings enter elsewhere in the
momentum update. In dry flow, :math:`q_{t,f}=0`.

Dividing :math:`B_z` by a representative density gives a buoyancy scale
with units of acceleration, but **the source alone is not a complete
prediction of vertical acceleration**.

Type 1: total-density perturbation
----------------------------------

The type-1 expression is the gravitational force associated with the
full total-density departure from the reference state:

.. math::

   B_{z,f}^{(1)}
   =g_z\,\overline{\left[\rho_d(1+q_t)-\rho_0(1+q_{v0})\right]}.

The entire bracketed density perturbation is averaged to the face. A
reference-state column with :math:`\rho_d=\rho_0` and
:math:`q_t=q_{v0}` has zero buoyancy. At otherwise fixed dry density,
adding condensate makes the air heavier and gives a downward contribution.
The compressible EOS also allows pressure perturbations to modify the
density, and type 1 retains that contribution.

Types 2 and 3: temperature perturbation
---------------------------------------

These two input values use the same approximate expression. Define

.. math::

   \epsilon_v=\frac{R_v}{R_d}-1,\qquad
   \beta_T=\frac{T-T_0}{T_0}
            +\epsilon_v(q_v-q_{v0})-q_{\mathrm{cond}},

where :math:`R_d` and :math:`R_v` are the gas constants of dry air and
water vapor, :math:`T` is diagnosed from the compressible thermodynamic
state, and :math:`T_0` is diagnosed from :math:`p_0` and :math:`\theta_0`.
In the active buoyancy helper, that reference conversion uses the fixed
exponent :math:`R_d/C_{p,d}`, with
:math:`C_{p,d}=1004.5\,\mathrm{J\,kg^{-1}\,K^{-1}}`.

The interior-face force is

.. math::

   B_{z,f}^{(2)}=B_{z,f}^{(3)}
   =-\overline{\rho_0}\,g_z\,\overline{\beta_T}.

The three contributions to :math:`\beta_T` have familiar interpretations:
warmer air favors upward buoyancy, additional water vapor favors upward
buoyancy through its gas constant, and condensate or precipitation mass
favors downward buoyancy. In a dry simulation, the water terms are zero.
This temperature-based expression is a **dilute-moisture,
small-pressure-perturbation approximation** to the density force, not an
identity for arbitrary compressible states.

Type 4: potential-temperature perturbation
------------------------------------------

Type 4 uses the dry potential-temperature anomaly instead of the actual
temperature anomaly. Define

.. math::

   \beta_\theta=\frac{\theta_d-\theta_0}{\theta_0}
           +\epsilon_v(q_v-q_{v0})-q_{\mathrm{cond}},

and compute

.. math::

   B_{z,f}^{(4)}
   =-\overline{\rho_0}\,g_z\,\overline{\beta_\theta}.

In dry flow this reduces to the potential-temperature perturbation term.
A warm potential-temperature perturbation gives an upward contribution;
a negative perturbation gives a downward contribution. As in the
temperature formulation, moisture enhances buoyancy through vapor
anomalies and reduces it through condensate loading.

This is an approximation to the full compressible density force. In
particular, a pressure perturbation at fixed potential temperature has
no explicit contribution to :math:`\beta_\theta`.

Anelastic buoyancy
------------------

In anelastic mode, ERF fixes the prognostic dry density to its base-state
value and imposes the anelastic velocity constraint. The active buoyancy
calculation does not use the compressible numeric selector.

For **dry anelastic flow**, ERF first averages the potential temperatures
to the vertical face and then forms the fractional perturbation:

.. math::

   B_{z,f}^{A,\mathrm{dry}}
   =-\overline{\rho_0}\,g_z\,
     \frac{\overline{\theta_d}-\overline{\theta_0}}
          {\overline{\theta_0}}.

For **moist anelastic flow**, ERF first forms
:math:`\beta_\theta` at each cell center and then averages it:

.. math::

   B_{z,f}^{A,\mathrm{moist}}
   =-\overline{\rho_0}\,g_z\,\overline{\beta_\theta}.

The order of averaging is different. Consequently, the two formulas do
not have to coincide exactly on a vertically varying reference profile
even when the water-mixing-ratio terms are zero. The pressure-gradient
term in anelastic momentum is supplied by the projection, not by
recomputing the compressible EOS pressure perturbation.

Why the compressible formulas differ
------------------------------------

The compressible equation of state implemented in ERF is

.. math::

   p=P_{00}
     \left[\frac{R_d\rho_d\theta_d(1+\alpha q_v)}{P_{00}}\right]^\Gamma,
   \qquad \alpha=\frac{R_v}{R_d},

where :math:`P_{00}=10^5\,\mathrm{Pa}` is the **fixed pressure used to
define potential temperature**, not the hydrostatic profile
:math:`p_0(z)`. The EOS exponent is fixed at :math:`\Gamma=1.4`. For
states consistent with this EOS, the exact total-density ratio is

.. math::

   \frac{\rho_t}{\rho_{t0}}
   =\left(\frac{p}{p_0}\right)^{1/\Gamma}
    \frac{\theta_0}{\theta_d}
    \frac{1+q_t}{1+q_{v0}}
    \frac{1+\alpha q_{v0}}{1+\alpha q_v}.

Writing :math:`p'=p-p_0`, :math:`\theta_d'=\theta_d-\theta_0`,
:math:`\delta q_v=q_v-q_{v0}`, and
:math:`\delta\rho_t=\rho_t-\rho_{t0}`, a first-order expansion gives

.. math::

   \frac{\delta\rho_t}{\rho_{t0}}
   =\frac{1}{\Gamma}\frac{p'}{p_0}
    -\frac{\theta_d'}{\theta_0}
    +\left[\frac{1}{1+q_{v0}}
           -\frac{\alpha}{1+\alpha q_{v0}}\right]\delta q_v
    +\frac{q_{\mathrm{cond}}}{1+q_{v0}}+O(\delta^2),

where :math:`O(\delta^2)` collects products and higher orders of small
state departures. Ignoring the pressure-density term and then taking the
dilute-water limit yields the type-4 expression. The currently
implemented approximation is **not** the exact finite-humidity
linearization, whose vapor and condensate coefficients depend on
:math:`q_{v0}`.

For small perturbations, with :math:`T'=T-T_0`,
:math:`T'/T_0\simeq\theta_d'/\theta_0+
(R_d/C_{p,d})p'/p_0`. Therefore types 2 and 3 retain a different
pressure-related response than type 4; neither is identical to the
full-density type 1.

.. note::

   A dry **pressure-only** perturbation illustrates the distinction.
   If :math:`p>p_0` while :math:`\theta_d=\theta_0`, type 1 gives a
   downward density-perturbation force, types 2/3 an upward
   temperature-perturbation force, and type 4 zero. These are separate
   **gravitational-source** contributions. They must not be interpreted
   as predictions of the complete vertical momentum tendency without
   also including the pressure-gradient term.

The runtime input ``erf.c_p`` defaults to the fixed reference heat
capacity :math:`C_{p,d}` but is a separate configurable parameter.
Changing ``erf.c_p`` affects selected temperature and forcing calculations but
**does not change** the fixed EOS :math:`\Gamma`. Consequently, results
that mix the runtime :math:`R_d/c_p` with the fixed EOS need not satisfy
all inverse thermodynamic identities away from the reference value.
See :ref:`GoverningEquations` and :ref:`ConstantsAndUnits`.

Supported configurations and further reading
--------------------------------------------

Embedded-boundary buoyancy has a separate implementation. With gravity
enabled, its established choices are **dry compressible type 1** and
**dry anelastic**; the other embedded-boundary combinations should not
be treated as supported buoyancy configurations. Compressible
``Kessler_NoRain``, ``SAM``, and ``SAM_NoPrecip_NoIce`` require type 1 in
the current source. These restrictions do not change the meaning of the
selector for other supported simulations.

For background on cloud-resolving models and the sensitivity of moist
models to their governing-equation choices, see
`Khairoutdinov and Randall (2003) <khairoutdinov2003cloud_>`_,
*Journal of the Atmospheric Sciences*, **60**, 607--625, and
`Bryan and Fritsch (2002) <bryan2002benchmark_>`_,
*Monthly Weather Review*, **130**, 2917--2928. These are background
references; the expressions above describe **ERF's implemented options**.

.. _khairoutdinov2003cloud: https://journals.ametsoc.org/view/journals/atsc/60/4/1520-0469_2003_060_0607_crmota_2.0.co_2.xml
.. _bryan2002benchmark: https://journals.ametsoc.org/view/journals/mwre/130/12/1520-0493_2002_130_2917_absfmn_2.0.co_2.xml
