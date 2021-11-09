
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _Algorithms:


Time Advance
============

Runge-Kutta
-----------

ERF uses the three-stage low-storage TVD RK3 scheme of `Gottlieb and Shu`_ to advance the solution.

.. _`Gottlieb and Shu`: https://www.ams.org/journals/mcom/1998-67-221/S0025-5718-98-00913-2/S0025-5718-98-00913-2.pdf

Specifically, for

.. math::

  \frac{d \mathbf{S}}{dt} = f(\mathbf{S})

where :math:`\mathbf{S}` is the solution vector, we solve

.. math::

  \mathbf{S}^{(1)} &=& \mathbf{S}^n + \Delta t f(\mathbf{S}^n)

  \mathbf{S}^{(2)} &=& \frac{3}{4} \mathbf{S}^n + \frac{1}{4} ( \mathbf{S}^{(1)} + \Delta t f(\mathbf{S}^{(1)}) )

  \mathbf{S}^{n+1} &=& \frac{1}{3} \mathbf{S}^n + \frac{2}{3} ( \mathbf{S}^{(2)} + \Delta t f(\mathbf{S}^{(2)}) )

In the code, the time-stepping is implemented in :cpp:`Source/RK3/RK3_driver.cpp`, while
:math:`f` is computed in :cpp:`Source/RK3/RK3_stage.cpp`

We note this can also be written as

.. math::

  \mathbf{S}^{(1)} &=& \mathbf{S}^n + \Delta t f(\mathbf{S}^n)

  \mathbf{S}^{(2)} &=& \mathbf{S}^n + \Delta t \; ( \frac{1}{4} f(\mathbf{S}^n) +  \frac{1}{4} f(\mathbf{S}^{(1)}) )

  \mathbf{S}^{n+1} &=& \mathbf{S}^n + \Delta t \; ( \frac{2}{3} f(\mathbf{S}^{(2)}) + \frac{1}{6} f(\mathbf{S}^{(1)}) +  \frac{1}{6} f(\mathbf{S}^{n}) )

which makes it clear that the boundary conditions for :math:`\mathbf{S}^{(1)}` live at time :math:`t^{n+1}`
while the boundary conditions for :math:`\mathbf{S}^{(2)}` are at time :math:`t^{n+1/2}`.

Acoustic Sub-stepping
---------------------

We sub-step the acoustic modes within each Runge-Kutta stage, following the algorithm
as described in `Klemp, Skamarock and Dudhia (2006)`_ but with further assumptions of height-based coordinates and a dry atmosphere.

.. _`Klemp, Skamarock and Dudhia (2006)`: https://journals.ametsoc.org/view/journals/mwre/135/8/mwr3440.1.xml

We first recall the equations in the form, here defining :math:`\mathbf{R}` for each equation to include all terms that contribute to the time evolution.

.. math::

  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}) + R_\rho,

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + p^\prime I) + {\mathbf F}_\mathbf{u} + \mathbf{R}_\mathbf{u}

  \frac{\partial (\Theta)}{\partial t} &=& - \nabla \cdot (\mathbf{u} \Theta) + \nabla \cdot (\alpha_{T}\ \nabla (\Theta)) + R_{\Theta},

  \frac{\partial (\rho C)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} C) + \nabla \cdot (\alpha_{C}\ \nabla (\rho C)) + R_{\rho C},

where we have defined :math:`\mathbf{U} = (U,V,W) = \rho \mathbf{u} = (\rho u, \rho v, \rho w)` , :math:`\Theta = \rho \theta` and
:math:`\mathbf{F}_\mathbf{U} = (F_U, F_V, F_W) = \rho^\prime \mathbf{g} + \nabla \cdot \tau + \mathbf{F}`

Using the relation :math:`\nabla p = \gamma R_d \pi \nabla \Theta,` where the Exner function :math:`\pi = (p/p_0)^\kappa` with :math:`\kappa = R_d / c_p,`
we can re-write the unapproximated momentum equations

.. math::

  \frac{\partial U}{\partial t} &=& - \nabla \cdot (\mathbf{u} U) - \gamma R_d \pi \frac{\partial \Theta^\prime}{\partial x} + F_u

  \frac{\partial V}{\partial t} &=& - \nabla \cdot (\mathbf{u} V) - \gamma R_d \pi \frac{\partial \Theta^\prime}{\partial y} + F_v

  \frac{\partial W}{\partial t} &=& - \nabla \cdot (\mathbf{u} W) - \gamma R_d \pi \frac{\partial \Theta^\prime}{\partial z}
                                                                              - g (\overline{\rho} \frac{\pi^\prime}{\overline{\pi}} - \rho^\prime) + F_w


We then define new perturbational quantities, e.g., :math:`\mathbf{U}^{\prime \prime} = \mathbf{U} - \mathbf{U}^t` where :math:`\mathbf{U}^t`
is the momentum at the most recent time reached by the "large" time step, which could be :math:`t^{n}` or one of the intermediate
Runge-Kutta steps.  Then the fully explicit acoustic substepping evolves the equations in the form

.. math::

  U^{\prime \prime, \tau + \delta \tau} - U^{\prime \prime, \tau} &=&  \delta \tau (
              -\gamma R_d \pi^t \frac{\partial \Theta^{\prime \prime, \tau}}{\partial x} + R^t_U)

  V^{\prime \prime, \tau + \delta \tau} - V^{\prime \prime, \tau} &=&  \delta \tau (
              -\gamma R_d \pi^t \frac{\partial \Theta^{\prime \prime, \tau}}{\partial y} + R^t_V)

  W^{\prime \prime, \tau + \delta \tau} - W^{\prime \prime, \tau} &=&  \delta \tau (
            -\gamma R_d \pi^t \frac{\partial \Theta^{\prime \prime, \tau}}{\partial z}
            - g \overline{\rho} \frac{R_d}{c_v} \frac{\pi^t}{\overline{\pi}} \frac{\Theta^{\prime \prime, t}}{\Theta^t}
            + g \rho^{\prime \prime, t} + R^t_W )

  \Theta^{\prime \prime, \tau + \delta \tau} - \Theta^{\prime \prime, \tau} &=& \delta \tau (
          -\frac{\partial (U^{\prime \prime, \tau} \theta^t)}{\partial x} +
          -\frac{\partial (V^{\prime \prime, \tau} \theta^t)}{\partial y} +
          -\frac{\partial (W^{\prime \prime, \tau} \theta^t)}{\partial z} +  R^t_{\Theta} )

  \rho^{\prime \prime, \tau + \delta \tau} -\rho^{\prime \prime, \tau} &=& \delta \tau (
          - \frac{\partial U^{\prime \prime, \tau}}{\partial x} - \frac{\partial V^{\prime \prime, \tau}}{\partial y}
          - \frac{\partial W^{\prime \prime, \tau}}{\partial z} +  R^t_{\rho} )

In the equations above the order of solution doesn't matter since the updates are fully explicit.

If we instead use the semi-implicit form below, then we must first update the horizontal momentum, then update the remaining coupled equations,
which requires a tri-diagonal solve.  The semi-implicit form of the equations have the same horizontal momentum equations but the
following for :math:`W^{\prime \prime}, \Theta^{\prime \prime}` and :math:`\rho^{\prime \prime}.`
We note that :math:`\beta = 1` below would correspond to fully explicit; :math:`\beta = 0` would be fully implicit for :math:`W^{\prime \prime}`.

.. math::

  W^{\prime \prime, \tau + \delta \tau} - W^{\prime \prime, \tau} &=&  \delta \tau (
            -\gamma R_d \pi^t \frac{\partial ( \beta \Theta^{\prime \prime, \tau}  (1 - \beta) \Theta^{\prime \prime, \tau + \delta \tau} ) }{\partial z} \\
            && - g \overline{\rho} \frac{R_d}{c_v} \frac{\pi^t}{\overline{\pi}}
             \frac{ ( \beta \Theta^{\prime \prime, \tau}  (1 - \beta) \Theta^{\prime \prime, \tau + \delta \tau} )}{\Theta^t}
            + g (\beta \rho^{\prime \prime, \tau} + (1 - \beta) \rho^{\prime \prime, \tau + \delta \tau } ) + R^t_W )

.. math::

  \Theta^{\prime \prime, \tau + \delta \tau} - \Theta^{\prime \prime, \tau} =  \delta \tau (
          -\frac{\partial (U^{\prime \prime, \tau + \delta \tau} \theta^t)}{\partial x}
          -\frac{\partial (V^{\prime \prime, \tau + \delta \tau} \theta^t)}{\partial y}
          -\frac{\partial (( \beta W^{\prime \prime, \tau} + (1 - \beta) W^{\prime \prime, \tau + \delta \tau} ) \theta^t)}{\partial z} +  R^t_{\Theta} )

.. math::

  \rho^{\prime \prime, \tau + \delta \tau} - \rho^{\prime \prime, \tau} =  \delta \tau (
          - \frac{\partial U^{\prime \prime, \tau + \delta \tau }}{\partial x}
          - \frac{\partial V^{\prime \prime, \tau + \delta \tau }}{\partial y}
          - \frac{\partial (\beta W^{\prime \prime, \tau} + (1-\beta) W^{\prime \prime, \tau + \delta \tau})}{\partial z} +  R^t_{\rho} )


We note that the only approximation in this system so far is in the linearization of the ideal gas law to define
:math:`\pi^{\prime \prime} = R_d \pi^t \Theta^{\prime \prime} / (c_v \Theta^t).`

Klemp et al note that with second-order differencing on a C grid, eliminating :math:`\rho^{\prime \prime}` and :math:`\Theta^{\prime \prime}`
from the vertical momentum equation using the final two equations results in a tridiagonal equation that is easily inverted.
