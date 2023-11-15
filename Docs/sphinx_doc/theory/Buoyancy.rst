
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Buoyancy:

Density of the mixture
========================

The total density in a given cell is given by

.. math::
    \rho &=& \frac{m}{V} = \frac{m_a + m_v + m_c + m_p}{V},

where :math:`m_a` is the mass of dry air, :math:`m_v` is the mass of water vapor, :math:`m_c` is the mass of liquid water, and :math:`m_p` is the mass of precipitate.
From the definitions of the mass mixing ratio (ratio of mass of a component to mass of dry air), we have for any component

.. math::
    q_i \equiv = \frac{m_i}{m_a}.

Using this we can write

.. math::
    \rho = m_a\frac{(1 + q_v + q_c + q_p)}{V}
          = \rho_d(1 + q_v + q_c + q_p),

where :math:`\rho_d \equiv \cfrac{m_a}{V}` is the density of dry air.


Buoyancy
=========

ERF has three options for how to define the buoyancy force.  Even in the absence of moisture these
expressions are not equivalent.

Type 1
------

One version of the buoyancy force is expressed simply as

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g}

.. math::
     \rho^\prime = \rho_{total} - \rho_0

where the full density :math:`\rho_{total}` is the sum of dry and moist components and :math:`\rho_0` is the base state density
for dry air only.

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
