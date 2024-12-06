
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _WetEquations:

Compressible Equations (Moist)
===============================

Model 1: Warm Moisture with no Precipitation
--------------------------------------------------

With this model, which is analogous to that in FASTEddy, we
consider a mixture of dry air :math:`\rho_d` and nonprecipitating water vapor :math:`\rho_v`,
assumed to be a perfect ideal gas with constant heat capacities
:math:`c_{vd}`, :math:`c_{vv}`, :math:`c_{pd}`, :math:`c_{pv}`, and
(non-precipitating) cloud water :math:`\rho_c`.

Neglecting the volume occupied by all water not in vapor form, we have

.. math::
  p = p_d + p_v = \rho_d R_d T + \rho_v R_v T

where :math:`p_d` and :math:`p_v` are the partial pressures of dry air and water vapor, respectively,
and :math:`R_d` and :math:`R_v` are the gas constants for dry air and water vapor, respectively.

We define the mixing ratio of each moist component, :math:`q_s`, as the mass density of species :math:`s`
relative to the density of dry air, i.e. :math:`q_s = \frac{\rho_s}{\rho_d}`.

Governing Equations
-------------------
The governing equations without precipitating moisture variables are

.. math::
   \frac{\partial \rho_d}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u}),

   \frac{\partial (\rho_d \mathbf{u})}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) - \frac{1}{1 + q_v + q_c} ( \nabla p^{\prime}  - \delta_{i,3}\mathbf{B} ) - \nabla \cdot \boldsymbol{\tau} + \mathbf{F}_{u},

   \frac{\partial (\rho_d \theta_d)}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \theta_d) + \nabla \cdot ( \rho_d \alpha_{\theta}\ \nabla \theta_d) + F_{\theta} + H_{n},

   \frac{\partial (\rho_d \boldsymbol{\phi})}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \boldsymbol{\phi}) + \nabla \cdot ( \rho_d \alpha_{\phi}\ \nabla \boldsymbol{\phi}) + \mathbf{F}_{\phi},

   \frac{\partial (\rho_d \mathbf{q_{n}})}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \mathbf{q_{n}}) + \nabla \cdot (\rho_d \alpha_{q} \nabla \mathbf{q_{n}}) + \mathbf{F_{n}},

the non-precipitating water mixing ratio vector :math:`\mathbf{q_{n}} = \left[ q_v \;\; q_c \;\; q_i \right]` includes water vapor, :math:`q_v`, cloud water, :math:`q_c`, and cloud ice, :math:`q_i`, although some models may not include cloud ice. The source terms for moisture variables, :math:`\mathbf{F_{n}}`, and their corresponding impact on potential temperature, :math:`H_{n}` are specific to the employed model. For the Kessler microphysics scheme, these terms are detailed in :ref:`sec:Kessler Microphysics model <Microphysics>`.


The pressure perturbation is computed as

.. math::
  p^\prime = P_{00} \left( \frac{R_d \rho_d \theta_m}{P_{00}} \right)^\gamma - p_{0}

where :math:`\gamma = c_{p} / (c_{p} - R_{d})` and

.. math::
  \theta_m = \theta_d (1 + \frac{R_v}{R_d} q_v)

is the moist potential temperature.  We note that this is the only place :math:`\theta_m` is used; we evolve :math:`\theta_d` above.

Model 2: Full Moisture Including Precipitation
--------------------------------------------------

With this model, in addition to dry air :math:`\rho_d` and nonprecipitating water vapor :math:`\rho_v`,
assumed to be a perfect ideal gas with constant heat capacities
:math:`c_{vd}`, :math:`c_{vv}`, :math:`c_{pd}`, :math:`C_{pv}`,
we include
non-precipitating condensates :math:`\rho_c + \rho_i`,
and precipitating condensates :math:`\rho_p = \rho_{rain} + \rho_{snow} + \rho_{graupel}`.
Here
:math:`\rho_c` is the density of cloud water and
:math:`\rho_i` is the density of cloud ice, and
we define the sum of all non-precipitating moist quantities to be :math:`\rho_T = \rho_v + \rho_c + \rho_i`.
All condensates  are treated as incompressible; cloud water and ice
have constant heat capacities :math:`C_p` and :math:`C_i`, respectively.

Neglecting the volume occupied by all water not in vapor form, we have

.. math::
  p = p_d + p_v = \rho_d R_d T + \rho_v R_v T

where :math:`p_d` and :math:`p_v` are the partial pressures of dry air and water vapor, respectively,
and :math:`R_d` and :math:`R_v` are the gas constants for dry air and water vapor, respectively.

We define the mixing ratio of each moist component, :math:`q_s`, as the mass density of species :math:`s`
relative to the density of dry air, i.e. :math:`q_s = \frac{\rho_s}{\rho_d}`.

We define the total potential temperature

.. math::
  \theta = \frac{\sum_s \rho_s \theta_s}{\sum_s \rho_s} \approx (\theta_d + q_v \theta_v + q_i \theta_i + q_c \theta_c).

and write the EOS as

.. math::
   T = \theta (\frac{p}{p_0})^\frac{R^\star}{C_p^\star}

or

.. math::
   p = P_{00} (\frac{\Pi}{c_p^\star})^{\frac{c_p^\star}{R^\star}}

where :math:`P_{00}` is the reference pressure. and

.. math::
  \Pi = C_p^\star (\frac{p}{\alpha p_0})^\frac{R^\star}{C_p^\star}

with :math:`\alpha = \frac{R^\star}{p}(\frac{p}{p_0})^\frac{R^\star}{c_p^\star} \theta`

here, :math:`R^\star =  R_{d} + q_v R_{v} + q_i R_{i} + q_p R_{p}`, and :math:`C_p^\star = C_{pd} + q_v C_{pv} + q_i C_{pi} + q_p C_{pp}`.

:math:`R_d`, :math:`R_v`, :math:`R_i`, and :math:`R_p` are the gas constants for dry air, water vapor, cloud ice, precipitating condensates, respectively. :math:`C_{pd}`, :math:`C_{pv}`, :math:`C_{pi}`, and :math:`C_{pp}` are the specific heats for dry air,
water vapor, cloud ice, and precipitating condensates, respectively.

Governing Equations
-------------------
The governing equations with precipitating moisture components are

.. math::
   \frac{\partial \rho_d}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u}),

   \frac{\partial (\rho_d \mathbf{u})}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) - \frac{1}{1 + q_v + q_c} ( \nabla p^{\prime}  - \delta_{i,3}\mathbf{B} ) - \nabla \cdot \boldsymbol{\tau} + \mathbf{F}_{u},

   \frac{\partial (\rho_d \theta_d)}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u}        \theta_d) + \nabla \cdot ( \rho_d \alpha_{\theta}\ \nabla \theta_d) + F_{\theta} + H_{n} + H_{p},

   \frac{\partial (\rho_d \boldsymbol{\phi})}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \boldsymbol{\phi}) + \nabla \cdot ( \rho_d \alpha_{\phi}\ \nabla \boldsymbol{\phi}) + \mathbf{F}_{\phi},

   \frac{\partial (\rho_d \mathbf{q_{n}})}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \mathbf{q_{n}}) + \nabla \cdot (\rho_d \alpha_{q} \nabla \mathbf{q_{n}}) + \mathbf{F_{n}} + \mathbf{G_{p}},

   \frac{\partial (\rho_d \mathbf{q_{p}})}{\partial t} = - \nabla \cdot (\rho_d \mathbf{u} \mathbf{q_{p}}) + \partial_{z} \left( \rho_d \mathbf{w_{t}} \mathbf{q_{p}} \right) + \mathbf{F_{p}}.

the non-precipitating water mixing ratio vector :math:`\mathbf{q_{n}} = \left[ q_v \;\; q_c \;\; q_i \right]` includes water vapor, :math:`q_v`, cloud water, :math:`q_c`, and cloud ice, :math:`q_i`, although some models may not include cloud ice; similarly, the precipitating water mixing ratio vector :math:`\mathbf{q_{p}} = \left[ q_r \;\; q_s \;\; q_g \right]` involves rain, :math:`q_r`, snow, :math:`q_s`, and graupel, :math:`q_g`, though some models may not include these terms. The source terms for moisture variables, :math:`\mathbf{F_{p}}`, :math:`\mathbf{F_{n}}`, :math:`\mathbf{G_{p}}`, and their corresponding impact on potential temperature, :math:`H_{n}` and :math:`H_{p}`, and the terminal velocity, :math:`\mathbf{w_{t}}` are specific to the employed model. For the Kessler microphysics scheme, these terms are detailed in :ref:`sec:Kessler Microphysics model <Microphysics>`.
