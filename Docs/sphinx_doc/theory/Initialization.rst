
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Initialization:

Initialization
==================

To initialize a background (base) state for a simulation with a hydrostatic atmosphere, the hydrostatic equation balancing the pressure gradient
and gravity must be satisfied. This section details the procedure for the initialization of the background state. The procedure is similar for the
cases with and without moisture, the only difference being that the density for the cases with moisture has to be the total density
:math:`\rho = \rho_d(1 + q_t)`, where :math:`\rho_d` is the dry density, and :math:`q_t` is the total mass mixing ratio -- water vapor and liquid water, instead
of the dry density :math:`\rho_d` for cases without moisture.

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

where :math:`\rho = \rho_d(1 + q_t)`, :math:`\rho_d` is the dry density, and :math:`q_t` is the total mass mixing ratio -- water vapor and liquid water. Using an average value of :math:`\rho` for the integration from :math:`z = z(k-1)` to :math:`z(k)`, we get

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
