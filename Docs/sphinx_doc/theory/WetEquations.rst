
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _WetEquations:

Prognostic Equations (Moist)
===============================

Model 1: Warm Moisture with no Precipitation
--------------------------------------------------

With this model, which is analogous to that in FASTEddy, we
consider a mixture of dry air :math:`\rho_d` and nonprecipitating water vapor :math:`\rho_v`,
assumed to be a perfect ideal gas with constant heat capacities
:math:`C_{vd}`, :math:`C_{vv}`, :math:`C_{pd}`, :math:`C_{pv}`, and
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
The governing equations for this model are

.. math::
  \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u})

  \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u}) -
          \frac{1}{1 + q_v + q_c} ( \nabla p^\prime  + \delta_{i,3}\mathbf{B} ) - \nabla \cdot \tau + \mathbf{F}

  \frac{\partial (\rho_d \theta_d)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \theta_d)
                + \nabla \cdot ( \rho_d \alpha_{T}\ \nabla \theta_d) + \frac{\theta_d L_v}{T_d C_{pd}} f_{cond}

  \frac{\partial (\rho_d q_v)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_vi) + \nabla \cdot (\rho_d \alpha \nabla q_v) - f_{cond}

  \frac{\partial (\rho_d q_c)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_c) + \nabla \cdot (\rho_d \alpha \nabla q_c) + f_{cond}

Here :math:`L_v` is the latent heat of vaporization, :math:`\theta_d` is the (dry) potential temperature
:math:`\mathbf{B}` is the buoyancy force, which is defined in :ref:`Buoyancy <Buoyancy>`.

The pressure perturbation is computed as

.. math::
  p^\prime = p_0 \left( \frac{R_d \rho_d \theta_m}{p_0} \right)^\gamma - p_0

where :math:`\gamma = C_{pd} / C_{vd}` and

.. math::
  \theta_m = \theta_d (1 + \frac{R_v}{R_d} q_v)

is the moist potential temperature.  We note that this is the only place :math:`\theta_m` is used; we evolve :math:`\theta_d` above.

Model 2: Full Moisture Including Precipitation
--------------------------------------------------

With this model, in addition to dry air :math:`\rho_d` and nonprecipitating water vapor :math:`\rho_v`,
assumed to be a perfect ideal gas with constant heat capacities
:math:`C_{vd}`, :math:`C_{vv}`, :math:`C_{pd}`, :math:`C_{pv}`,
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
   p = p_0 (\frac{\Pi}{C_p^\star})^{\frac{C_p^\star}{R^\star}}

where :math:`p_0` is the reference pressure. and

.. math::
  \Pi = C_p^\star (\frac{p}{\alpha p_0})^\frac{R^\star}{C_p^\star}

with :math:`\alpha = \frac{R^\star}{p}(\frac{p}{p_0})^\frac{R^\star}{c_p^\star} \theta`

here, :math:`R^\star =  R_{d} + q_v R_{v} + q_i R_{i} + q_p R_{p}`, and :math:`C_p^\star = C_{pd} + q_v C_{pv} + q_i C_{pi} + q_p C_{pp}`.

:math:`R_d`, :math:`R_v`, :math:`R_i`, and :math:`R_p` are the gas constants for dry air, water vapor, cloud ice, precipitating condensates, respectively. :math:`C_{pd}`, :math:`C_{pv}`, :math:`C_{pi}`, and :math:`C_{pp}` are the specific heats for dry air,
water vapor, cloud ice, and precipitating condensates, respectively.

Governing Equations
-------------------
We assume that all species have same average speed,
Then the governing equations become

.. math::
  \frac{\partial \rho_d}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} + \mathbf{F}_\rho)

  \frac{\partial (\rho_d \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \mathbf{u} + \mathbf{F}_u) -
          \frac{1}{1 + q_T + q_p}  \nabla p^\prime - \nabla \cdot \tau + \mathbf{F} + \delta_{i,3}\mathbf{B}

  \frac{\partial (\rho_d \theta)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} \theta + F_{\theta}) + \nabla \cdot ( \rho_d \alpha_{T}\ \nabla \theta) + F_Q

  \frac{\partial (\rho_d q_T)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_T +F_{q_{T}}) - Q

  \frac{\partial (\rho_d q_p)}{\partial t} &= - \nabla \cdot (\rho_d \mathbf{u} q_p + F_{q_{p}}) + Q

In this set of equations, the subgrid turbulent parameterization effects are included with fluxes
:math:`F_\rho`, :math:`F_u`, :math:`F_C`, :math:`F_{\theta}`, :math:`F_{q_{T}}`, :math:`F_{q_{r}}`.
:math:`\mathbf{F}` stands for the external force, and :math:`Q` and :math:`F_Q` represent the mass and energy transformation
of water vapor to/from water through condensation/evaporation, which is determined by the microphysics parameterization processes.
:math:`\mathbf{B}` is the buoyancy force, which is defined in :ref:`Buoyancy <Buoyancy>`.
