
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Buoyancy:

Buoyancy
=========

ERF has three options for how to define the buoyancy force.  Even in the absence of moisture these
expressions are not equivalent.

Density of the mixture
-----------------------

The total density in a cell containing air, water vapor, liquid water and precipitates is given by

.. math::
    \rho = \frac{m}{V} = \frac{m_a + m_v + m_c + m_p}{V},

where :math:`m_a` is the mass of dry air, :math:`m_v` is the mass of water vapor, :math:`m_c` is the mass of liquid water, and :math:`m_p` is the mass of precipitate.
From the definitions of the mass mixing ratio (ratio of mass of a component to mass of dry air), we have for any component

.. math::
    q_i \equiv \frac{m_i}{m_a}.

Using this we can write

.. math::
    \rho = m_a\frac{(1 + q_v + q_c + q_p)}{V}
          = \rho_d(1 + q_v + q_c + q_p),

where :math:`\rho_d \equiv \cfrac{m_a}{V}` is the density of dry air.


Type 1
------

One version of the buoyancy force is expressed simply as

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g}

.. math::
     \rho^\prime = \rho_{total} - \rho_0

where the total density :math:`\rho_{total} = \rho_d(1 + q_v + q_c + q_p)` is the sum of dry and moist components and :math:`\rho_0` is the total density
for the background state. For eg., a usual scenario is that of a background state that contains only air and vapor and no cloud water or precipitates. For such a state,
the total background density :math:`\rho_0 = \rho_{d_0}(1 + q_{v_0})`, where :math:`\rho_{d_0}` and :math:`q_{v_0}` are the background dry density and vapor mixing ratio respectively.
As a check, we observe that :math:`\rho^\prime_0 = 0`, which means that the background state is not buoyant.

Type 2
------

The second option for the buoyancy force is

.. math::
   \mathbf{B} = -\rho_0 \mathbf{g} ( 0.61 q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime
                  + \frac{T^\prime}{\bar{T}} (1.0 + 0.61 \bar{q_v} - \bar{q_i} - \bar{q_c} - \bar{q_p}) )

To derive this expression, we define :math:`T_v = T (1 + 0.61 q_v − q_c − q_i - q_p)`, then we can write

.. math::
    p = \rho (R_d q_d + R_v q_v) T = \rho R_d T (1 + 0.61 q_v − q_c − q_i - q_p ) = \rho R_d T_v


Starting from :math:`p = \rho R_d T_v` and neglecting :math:`\frac{p^\prime}{\bar{p}}`, we now write

.. math::
   \frac{\rho^\prime}{\overline{\rho}} = -\frac{T_v^\prime}{\overline{T_v}}

and define

.. math::

   T_v^\prime = T_v - \overline{T_v} \approx \overline{T} [ 0.61 q_v^\prime - (q_c^\prime + q_i^\prime + q_p^\prime)] +
               (T - \overline{T}) [1+ 0.61 \bar{q_v} - \bar{q_c} - \bar{q_i} - \bar{q_p} ] .

where we have retained only first order terms in perturbational quantities.

Then

.. math::

   \mathbf{B} = \rho^\prime \mathbf{g} = -\overline{\rho} \frac{\overline{T}}{\overline{T_v}} \mathbf{g} [ 0.61 q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime ) + \frac{T^\prime}{\overline{T_v}} (1.0 + 0.61 \bar{q_v} - \bar{q_i} - \bar{q_c} - \bar{q_p}) ]

where the overbar represents a horizontal average of the current state and the perturbation is defined relative to that average.

Again keeping only the first order terms in the mass mixing ratios, we can simplify this to

.. math::
   \mathbf{B} = \rho^\prime \mathbf{g} = -\rho_0 \mathbf{g} [ 0.61 q_v^\prime - q_c^\prime + q_i^\prime + q_p^\prime
                  + \frac{T^\prime}{\overline{T}} (1.0 + 0.61 \bar{q_v} - \bar{q_i} - \bar{q_c} - \bar{q_p}) ]

We note that this reduces to Type 3 if the horizontal averages of the moisture terms are all zero.

Type 3
------

The third formulation of the buoyancy term assumes that the horizontal averages of the moisture quantities are negligible,
which removes the need to compute horizontal averages of these quantities.   This reduces the Type 2 expression to the following:

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g} \approx -\rho_0 \mathbf{g} ( \frac{T^\prime}{\overline{T}}
                 + 0.61 q_v - q_c - q_i - q_p)

We note that this version of the buoyancy force matches that given in Marat F. Khairoutdinov and David A. Randall's paper (J. Atm Sciences, 607, 1983)
if we neglect :math:`\frac{p^\prime}{\bar{p_0}}`.

Type 4
------
.. math::

    \begin{equation}
    \mathbf{B} = \rho'\mathbf{g} \approx -\rho\Bigg(\frac{T'}{T} + 0.61 q_v' - q_c - q_p - \frac{p'}{p}\Bigg),
    \end{equation}

The derivation follows. The total density is given by :math:`\rho = \rho_d(1 + q_v + q_c + q_p)`, which can be written as

.. math::

    \rho = \frac{p (1 + q_v + q_c + q_p)}{R_dT\Bigg(1 + \cfrac{R_v}{R_d}q_v\Bigg)}

This can be written using binomial expansion as

.. math::

    \begin{align*}
    \rho &= \frac{p}{R_dT} (1 + q_v + q_c + q_p)\Bigg(1 + \frac{R_v}{R_d}q_v\Bigg)^{-1} \\
    &= \frac{p}{R_dT} (1 + q_v + q_c + q_p)\Bigg(1 - \frac{R_v}{R_d}q_v + O(q_v^2)\Bigg) \\
    &= \frac{p}{R_dT}\Bigg(1 + q_v + q_c + q_p - \frac{R_v}{R_d}q_v +  \text{H.O.T. such as } O(q_v^2) + O(q_vq_c)\Bigg) \\
    &\approx \frac{p}{R_dT}\Bigg(1 + q_v + q_c + q_p - \frac{R_v}{R_d}q_v\Bigg)
    \end{align*}

Taking log on both sides, we get

.. math::

    \log{\rho} = \log{p} - \log{R_d} - \log{T} + \log(1 - 0.61 q_v + q_c + q_p)

Taking derivative gives

.. math::

    \frac{\rho'}{\rho} = \frac{p'}{p} - \frac{T'}{T} + \frac{(-0.61 q_v' + q_c' + q_p')}{(1 - 0.61 q_v + q_c + q_p)}

Using :math:`- 0.61 q_v + q_c + q_p \ll 1`, we have

.. math::

    \frac{\rho'}{\rho} = \frac{p'}{p} - \frac{T'}{T} + (-0.61 q_v' + q_c' + q_p')

Since the background values of cloud water and precipitate mass mixing ratios -- :math:`q_c` and :math:`q_p` are zero, we have :math:`q_c' = q_c` and :math:`q_p' = q_p`. Hence, we have

.. math::

    \begin{equation}
    \rho'\approx -\rho\Bigg(\frac{T'}{T} + 0.61 q_v' - q_c - q_p - \frac{p'}{p}\Bigg),
    \end{equation}

