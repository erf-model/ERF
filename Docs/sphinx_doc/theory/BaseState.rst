
 .. role:: cpp(code)
    :language: c++

.. _sec:BaseState:

Construction of hydrostatically stratified base state
=========================================================

Here we describe how ERF initializes base state values of density
and potential temperature such that the density, pressure and potential
temperature satisfy both the hydrostatic balance and the equation of state.
The base-state density stored in ERF is the dry-air density
:math:`\rho_{d,0}`. The only moisture variable carried in the base state is
the water vapor mixing ratio :math:`q_{v,0}`; when total base-state density is
needed, ERF forms :math:`\rho_{d,0}(1 + q_{v,0})`.

Users have the option to define a dry or moist background state.

Computation of the dry density
-------------------------------
We express the total pressure :math:`p` as:

.. math::
   p = \rho_d R_d T_v.

By definition, we have:

.. math::
   T_v = \theta_m\left(\frac{p}{p_0}\right)^{R_d/C_p},

.. math::
   T = \theta_d\left(\frac{p}{p_0}\right)^{R_d/C_p}.

This gives:

.. math::
   p = \rho_d R_d \theta_m\left(\frac{p}{p_0}\right)^{R_d/C_p},

which, on rearranging, yields:

.. math::
   p = p_0\left(\frac{\rho_d R_d \theta_m}{p_0}\right)^\gamma.

To obtain :math:`\theta_m`, consider the density of dry air:

.. math::
   \rho_d = \frac{p - p_v}{R_d T}.

Substituting for :math:`\rho_d` from the above equation, we get:

.. math::
   \frac{p}{T_v} = \frac{p - p_v}{T},

which implies:

.. math::
   \frac{T_v}{T} = \frac{p}{p - p_v} = \frac{1}{1-\left(\cfrac{p_v}{p}\right)}.

We also have:

.. math::
   q_v = \frac{\rho_v}{\rho_d} = \frac{p_v}{R_v T}\frac{R_d T}{p-p_v} = \frac{r p_v}{p - p_v},

where :math:`p_v` is the partial pressure of water vapor, :math:`r = R_d/R_v \approx 0.622`, and :math:`q_v` is the vapor mass mixing ratio—the ratio of the
mass of vapor to the mass of dry air. Rearranging and using :math:`q_v \ll r`, we get:

.. math::
   \frac{p_v}{p} = \frac{1}{1 + \left(\cfrac{r}{q_v}\right)} \approx \frac{q_v}{r},

which, on substitution in the equation for :math:`\frac{T_v}{T}`, gives:

.. math::
   \frac{T_v}{T} = \frac{1}{1 - \left(\cfrac{q_v}{r}\right)}.

As :math:`q_v \ll 1`, a binomial expansion, ignoring higher-order terms, gives:

.. math::
   T_v = T\left(1 + \frac{R_v}{R_d}q_v\right).

Hence, the density of dry air is given by:

.. math::
   \rho_d = \frac{p}{R_d T_v} = \frac{p}{R_d T\left(1 + \cfrac{R_v}{R_d}q_v\right)}.


Initialization with a second-order integration of the hydrostatic equation
----------------------------------------------------------------------------

We have the hydrostatic equation given by

.. math::

    \frac{\partial p}{\partial z} = -\rho g,

where :math:`\rho = \rho_d(1 + q_v)`, :math:`\rho_d` is the dry density, and :math:`q_v` is the water vapor mixing ratio carried in the base state. Using an average value of :math:`\rho` for the integration from :math:`z = z(k-1)` to :math:`z(k)`, we get

.. math::

    p(k) = p(k-1) - \frac{(\rho(k-1) + \rho(k))}{2} g\Delta z.

The density at a point is a function of the pressure, potential temperature, and relative humidity. The latter two quantities are computed using user-specified profiles, and hence, for simplicity, we write :math:`\rho(k) = f(p(k))`. Hence

.. math::

    p(k) = p(k-1) - \frac{\rho(k-1)}{2}g\Delta z - \frac{f(p(k))}{2}g\Delta z.

Now, we define

.. math::

    F(p(k)) \equiv p(k) - p(k-1) + \frac{\rho(k-1)}{2}g\Delta z + \frac{f(p(k))}{2}g\Delta z = 0.

This is a non-linear equation in :math:`p(k)`. Consider a Newton-Raphson iteration (where :math:`n` denotes the iteration number) procedure

.. math::

    F(p+\delta p) \approx F(p) + \delta p \frac{\partial F}{\partial p} = 0,

which implies

.. math::

    \delta p = -\frac{F}{F'},

with the gradient being evaluated as

.. math::

    F' = \frac{F(p+\epsilon) - F(p)}{\epsilon},

and the iteration update is given by

.. math::

    p^{n+1} = p^n + \delta p.

For the first cell (:math:`k=0`), which is at a height of :math:`z = \frac{\Delta z}{2}` from the base, we have

.. math::

    p(0) = p_0 - \rho(0)g\frac{\Delta z}{2},

where :math:`p_0 = 1e5 \, \text{N/m}^2` is the pressure at the base. Hence, we define

.. math::

    F(p(0)) \equiv p(0) - p_0 + \rho(0)g\frac{\Delta z}{2},

and the Newton-Raphson procedure is the same.

.. _subsec:base-state-refined:

The base state on a refined mesh
----------------------------------------------------------------------------

With mesh refinement, the base state is constructed **independently on each level**,
against that level's own cell-centered heights :math:`z_{cc}`, by the same procedure
described above.  Running the same construction on every level, rather than
interpolating the coarse answer, is what leaves each level in discrete hydrostatic
balance on its own mesh.

Before that construction runs, the base state on a refined level is filled by
conservative interpolation from its parent, and the physical boundary conditions for the
base state are applied.  That is what gives values to the part of the refined level that
lies inside the domain but outside the refined grids, and to the ghost cells; the
per-level construction then overwrites only the refined grids.

Two consequences are worth stating explicitly.

**A refined level need not span the depth of the domain.**  The vertical integration is
causal upward: the value in a cell depends only on the cells at or below it.  A refined
region that stops partway up the domain therefore gets exactly the base state it would
have had if it had been refined all the way to the model top -- the values in the cells it
does contain are unchanged.  What the integration does require is a valid starting value
in the lowest cell of each column, which is why the grids must not be decomposed in the
vertical (see :ref:`subsec:no-vertical-decomposition`); ``amr.refine_grid_layout_z``
defaults to 0 in ERF, so this holds unless it is overridden.

**The base state does not depend on the coupling type.**  Under two-way coupling ERF
averages :math:`\det J` down from fine to coarse, so that the volume weighting used when
averaging the state down telescopes and remains conservative.  The cell-centered heights
:math:`z_{cc}` are deliberately **not** averaged down.  Were they, the coarse heights
would change after the coarse base state had already been built against them, leaving that
base state out of hydrostatic balance with the heights the dynamics subsequently use for
vertical gradients, Rayleigh damping and the sponge zones.  Not averaging them also keeps
:math:`z_{cc}` consistent with the nodal heights :math:`z_{nd}` it is derived from, which
are not averaged down either.

Note that the base states of two levels will not be identical to one another wherever
their terrain differs -- for instance when a finer level resolves topography that its
parent does not, either because it re-reads a terrain file at its own resolution or
because it reads a nested wrfinput file.  Each level is individually hydrostatic; the
difference between them is set by the difference in surface elevation and is carried up
the column as an essentially constant offset in pressure.
