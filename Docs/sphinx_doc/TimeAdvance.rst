
 .. role:: cpp(code)
    :language: c++

 .. _TimeAdvance:

Time Advance
============

Compressible Advance
---------------------

To advance the fully compressible solution in time, ERF uses a 3rd order Runge-Kutta method.
By default, acoustic substepping which solves the perturbational equations
implicitly in the vertical direction is used in each Runge-Kutta stage,
following the approach of `Klemp, Skamarock and Dudhia (2007)`_
However, there is a run-time option to turn off the substepping completely,
or to use an explicit rather than implicit solve in the substepping.

.. _`Klemp, Skamarock and Dudhia (2007)`: https://journals.ametsoc.org/view/journals/mwre/135/8/mwr3440.1.xml

Specifically, the 3rd order Runge-Kutta method solves

.. math::

  \frac{d \mathbf{S}}{dt} = f(\mathbf{S})

where :math:`\mathbf{S}` is the solution vector, in the following three steps:

.. math::

  \mathbf{S}^{*}   &=& \mathbf{S}^n + \frac{1}{3} \Delta t f(\mathbf{S}^n)

  \mathbf{S}^{**}  &=& \mathbf{S}^n + \frac{1}{2} \Delta t f(\mathbf{S}^{*})

  \mathbf{S}^{n+1} &=& \mathbf{S}^n +             \Delta t f(\mathbf{S}^{**})

.. _AnelasticTimeAdvance:

Anelastic Advance
---------------------

When solving the anelastic rather than fully compressible equations, ERF uses a 2nd order Runge-Kutta method
(with no substepping):

Specifically, the 2nd order Runge-Kutta method solves

.. math::

  \frac{d \mathbf{S}}{dt} = f(\mathbf{S})

where :math:`\mathbf{S}` is the solution vector, in the following two steps:

.. math::

  \mathbf{S}^{*}   &=& \mathbf{S}^n + \Delta t f(\mathbf{S}^n)

  \mathbf{S}^{n+1} &=& \mathbf{S}^n + \frac{\Delta t}{2} ( f(\mathbf{S}^{n}) + f(\mathbf{S}^{*}) )

.. _AcousticSubstep:

Acoustic Sub-stepping
---------------------

When solving the fully compressible equation set, by default we substep the acoustic modes within each Runge-Kutta stage.

Recall the equations in the following form,
here defining :math:`\mathbf{R}` for each equation to include all additional terms that contribute to the time evolution.

.. math::

  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}) = R_\rho,

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + p^\prime I) + {\mathbf F}_\mathbf{u} = \mathbf{R}_\mathbf{u}

  \frac{\partial (\Theta)}{\partial t} &=& - \nabla \cdot (\mathbf{u} \Theta) + \rho \alpha_{T}\ \nabla^2 \theta = R_{\Theta},

  \frac{\partial (\rho C)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} C) + \rho \alpha_{C}\ \nabla^2 C = R_{\rho C},

where we have defined :math:`\mathbf{U} = (U,V,W) = \rho \mathbf{u} = (\rho u, \rho v, \rho w)` , :math:`\Theta = \rho \theta` and
:math:`\mathbf{F}_\mathbf{U} = (F_U, F_V, F_W) = \rho^\prime \mathbf{g} + \nabla \cdot \tau + \mathbf{F}`

Using the relation :math:`\nabla p = \gamma R_d \pi \nabla \Theta,` where the Exner function :math:`\pi = (p/p_0)^\kappa` with :math:`\kappa = R_d / c_p,`
we can re-write the unapproximated momentum equations

.. math::

  \frac{\partial U}{\partial t} &=& - \nabla \cdot (\mathbf{u} U) - \gamma R_d \pi \frac{\partial \Theta^\prime}{\partial x} + F_u

  \frac{\partial V}{\partial t} &=& - \nabla \cdot (\mathbf{u} V) - \gamma R_d \pi \frac{\partial \Theta^\prime}{\partial y} + F_v

  \frac{\partial W}{\partial t} &=& - \nabla \cdot (\mathbf{u} W) - \gamma R_d \pi \frac{\partial \Theta^\prime}{\partial z}
                                                                              - g (\overline{\rho} \frac{\pi^\prime}{\overline{\pi}} - \rho^\prime) + F_w


We then define new perturbational quantities, e.g., :math:`\mathbf{U}^{\prime \prime} = \mathbf{U} - \mathbf{U}^t`
where :math:`\mathbf{U}^t` is the momentum at the most recent time reached by the "large" time step,
which could be :math:`t^{n}` or one of the intermediate Runge-Kutta steps.
Then the acoustic substepping evolves the equations in the form

.. math::

  U^{\prime \prime, \tau + \delta \tau} - U^{\prime \prime, \tau} &= \delta \tau \left(
              -\gamma R_d \pi^t \frac{\partial \Theta^{\prime \prime, \tau}}{\partial x} + R^t_U
              \right)

  V^{\prime \prime, \tau + \delta \tau} - V^{\prime \prime, \tau} &= \delta \tau \left(
              -\gamma R_d \pi^t \frac{\partial \Theta^{\prime \prime, \tau}}{\partial y} + R^t_V
              \right)

.. math::

  W^{\prime \prime, \tau + \delta \tau} - W^{\prime \prime, \tau} = \delta \tau \biggl(
          &-\gamma R_d \pi^t \frac{\partial (\beta_1 \Theta^{\prime \prime, \tau} +
                                              \beta_2 \Theta^{\prime \prime, \tau  + \delta \tau} ) }{\partial z} \\
          & - g \overline{\rho} \frac{R_d}{c_v} \frac{\pi^t}{\overline{\pi}}
             \frac{ (\beta_1 \Theta^{\prime \prime, \tau}  +
                     \beta_2 \Theta^{\prime \prime, \tau + \delta \tau} )}{\Theta^t} \\
          & + g (\beta_1 \rho^{\prime \prime, \tau} + \beta_2 \rho^{\prime \prime, \tau + \delta \tau } ) \\
          & + R^t_W \biggr)

.. math::

  \Theta^{\prime \prime, \tau + \delta \tau} - \Theta^{\prime \prime, \tau} =  \delta \tau \left(
          -\frac{\partial (U^{\prime \prime, \tau + \delta \tau} \theta^t)}{\partial x}
          -\frac{\partial (V^{\prime \prime, \tau + \delta \tau} \theta^t)}{\partial y}
          -\frac{\partial \left(( \beta_1 W^{\prime \prime, \tau} + \beta_2 W^{\prime \prime, \tau + \delta \tau} ) \theta^t\right)}{\partial z} +  R^t_{\Theta}
          \right)

.. math::

  \rho^{\prime \prime, \tau + \delta \tau} - \rho^{\prime \prime, \tau} =  \delta \tau \left(
          - \frac{\partial U^{\prime \prime, \tau + \delta \tau }}{\partial x}
          - \frac{\partial V^{\prime \prime, \tau + \delta \tau }}{\partial y}
          - \frac{\partial (\beta_1 W^{\prime \prime, \tau} + \beta_2 W^{\prime \prime, \tau + \delta \tau})}{\partial z} +  R^t_{\rho}
            \right)

where :math:`\beta_1 = 0.5 (1 - \beta_s)` and :math:`\beta_2 = 0.5 (1 + \beta_s)` with :math:`\beta_s = 0.1`.
:math:`\beta_s` is the acoustic step off-centering coefficient and 0.1 is the typical WRF value. This off-centering is intended to provide damping of both horizontally and vertically propagating sound waves by biasing the time average toward the future time step.

To solve the coupled system, we first evolve the equations for :math:`U^{\prime \prime, \tau + \delta \tau}`  and
:math:`V^{\prime \prime, \tau + \delta \tau}` explicitly using :math:`\Theta^{\prime \prime, \tau}` which is already known.
We then solve a tridiagonal system for :math:`W^{\prime \prime, \tau + \delta \tau}`, and once :math:`W^{\prime \prime, \tau + \delta \tau}`
is known, we update :math:`\rho^{\prime \prime, \tau + \delta \tau}` and :math:`\Theta^{\prime \prime, \tau + \delta \tau}.`

In addition to the acoustic off-centering, divergence damping is also applied
to control horizontally propagating sound waves.

.. math::

   p^{\prime\prime,\tau*} = p^{\prime\prime,\tau}
     + \beta_d \left( p^{\prime\prime,\tau} + p^{\prime\prime,\tau-\delta\tau} \right)

where :math:`\tau*` is the forward projected value used in RHS of the acoustic
substepping equations for horizontal momentum. According to Skamarock et al,
This is equivalent to including a horizontal diffusion term in the continuity
equation. A typical damping coefficient of :math:`\beta_d = 0.1` is used, as in
WRF.
