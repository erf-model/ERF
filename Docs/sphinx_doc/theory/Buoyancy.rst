.. role:: cpp(code)
   :language: c++

.. role:: f(code)
   :language: fortran

.. _Buoyancy:

Buoyancy
=========

ERF has several options for how to define the buoyancy force. All models may be employed
with the compressible formulation but those applicable to the anelastic formulation will be
explicitly stated in the description.

The buoyancy formulation is selected via the ``erf.buoyancy_type`` parameter. See :ref:`sec:Inputs`
for available values and defaults.

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

Type 1: Density Perturbation
-----------------------------

One version of the buoyancy force is expressed simply as

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g}

.. math::
     \rho^\prime = \rho_{total} - \rho_0

where the total density :math:`\rho_{total} = \rho_d(1 + q_v + q_c + q_p)` is the sum of dry and moist components and :math:`\rho_0` is the total density
for the background state. For example, a usual scenario is that of a background state that contains only air and vapor and no cloud water or precipitates. For such a state,
the total background density :math:`\rho_0 = \rho_{d_0}(1 + q_{v_0})`, where :math:`\rho_{d_0}` and :math:`q_{v_0}` are the background dry density and vapor mixing ratio respectively.
As a check, we observe that :math:`\rho^\prime_0 = 0`, which means that the background state is not buoyant.

This is the default formulation for both dry and certain moist simulations.

Type 2/3: Temperature Perturbation
-----------------------------------

**Note:** Types 2 and 3 are implemented identically in the code.

For **dry** simulations, the buoyancy is:

.. math::
     \mathbf{B} = -\rho_0 \mathbf{g} \frac{T'}{T_0}

For **moist** simulations, this formulation assumes that the horizontal averages of the moisture quantities are negligible:

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g} \approx -\rho_0 \mathbf{g} \left( \frac{T^\prime}{\overline{T}}
                 + 0.61 q_v - q_c - q_i - q_p \right)

We note that this version of the buoyancy force matches that given in Marat F. Khairoutdinov and David A. Randall's paper (J. Atm Sciences, 607, 1983)
if we neglect :math:`\frac{p^\prime}{\bar{p_0}}`.

Type 3 is utilized when the anelastic formulation is employed. A specialized implementation
based on dry potential temperature perturbations is used in the anelastic case.

Type 4: Potential Temperature Perturbation
-------------------------------------------

This expression for buoyancy is from `khairoutdinov2003cloud`_ and `bryan2002benchmark`_.

.. _`khairoutdinov2003cloud`: https://journals.ametsoc.org/view/journals/atsc/60/4/1520-0469_2003_060_0607_crmota_2.0.co_2.xml
.. _`bryan2002benchmark`: https://journals.ametsoc.org/view/journals/mwre/130/12/1520-0493_2002_130_2917_absfmn_2.0.co_2.xml

For **dry** simulations:

.. math::
    \mathbf{B} = -\rho_0 \mathbf{g} \frac{\theta'}{\theta_0}

    \begin{equation}
    \mathbf{B} = \rho'\mathbf{g} \approx -\rho \mathbf{g} \Bigg(\frac{T'}{T} + 0.61 q_v' - q_c - q_p - \frac{p'}{p}\Bigg)
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

which gives the final expression for Type 4 moist buoyancy.

Type 5: Anelastic Formulation (Internal Use Only)
--------------------------------------------------

.. note::
   Type 5 is not user-selectable via ``erf.buoyancy_type``. It describes the specialized
   buoyancy implementation used internally when the anelastic formulation is enabled.

Utilizing :math:`\theta_d` and neglecting the pressure term in Type 4 leads to:

.. math::

    \begin{equation}
    \rho'\approx -\rho\left(\frac{\theta_d'}{\theta} + 0.61 q_v' - q_c' - q_p'\right).
    \end{equation}

This buoyancy model is employed when utilizing the anelastic formulation.
