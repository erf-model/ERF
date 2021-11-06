
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

  \frac{d \mathbf{U}}{dt} = f(\mathbf{U})

where :math:`\mathbf{U}` is the solution vector, we solve

.. math::

  \mathbf{U}^{(1)} = \mathbf{U}^n + \Delta t f(\mathbf{U}^n)

  \mathbf{U}^{(2)} = \frac{3}{4} \mathbf{U}^n + \frac{1}{4} ( \mathbf{U}^{(1)} + \Delta t f(\mathbf{U}^{(1)}) )

  \mathbf{U}^{n+1} = \frac{1}{3} \mathbf{U}^n + \frac{2}{3} ( \mathbf{U}^{(2)} + \Delta t f(\mathbf{U}^{(2)}) )

In the code, the time-stepping is implemented in :cpp:`Source/RK3/RK3_driver.cpp`, while
:math:`f` is computed in :cpp:`Source/RK3/RK3_stage.cpp`

We note this can also be written as

.. math::

  \mathbf{U}^{(1)} = \mathbf{U}^n + \Delta t f(\mathbf{U}^n)

  \mathbf{U}^{(2)} = \mathbf{U}^n + \Delta t ( \frac{1}{4} f(\mathbf{U}^n) +  \frac{1}{4} f(\mathbf{U}^{(1)}) )

  \mathbf{U}^{n+1} = \mathbf{U}^n + \Delta t ( \frac{2}{3} \mathbf{U}^{(2)}  +  \frac{1}{6} \mathbf{U}^{(1)} +  \frac{1}{6} \mathbf{U}^{n} )

which makes it clear that the boundary conditions for :math:`\mathbf{U}^{(1)}` live at time :math:`t^{n+1}`
while the boundary conditions for :math:`\mathbf{U}^{(2)}` are at time :math:`t^{n+1/2}`.

Acoustic Sub-stepping
---------------------

In order to use the advective rather than acoustic CFL constraint to determine the time step, we sub-step
the acoustic modes within each Runge-Kutta stage, holding the non-acoustic terms fixed.  We follow the algorithm
as described in `Klemp, Skamarock and Dudhia (2006)`_

.._`Klemp, Skamarock and Dudhia (2006)`: https://journals.ametsoc.org/view/journals/mwre/135/8/mwr3440.1.xml

We first recall the equations in the form, here defining $\mathbf{R}$ for each equation to include all terms that contribute to the time evolution.

.. math::

  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}) = R_\rho,

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + p^\prime I) + {\mathbf F}_\mathbf{u} = \mathbf{R}_\mathbf{u}

  \frac{\partial (\Theta)}{\partial t} &=& - \nabla \cdot (\mathbf{u} \Theta) + \nabla \cdot (\alpha_{T}\ \nabla (\Theta)) = R_{\Theta},

  \frac{\partial (\rho C)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} C) + \nabla \cdot (\alpha_{C}\ \nabla (\rho C)) = R_{\rho C},

where we have defined :math:`\Theta = \rho \theta,` and
:math:`\mathbf{F}_\mathbf{u} = (F_u, F_v, F_w) = \rho^\prime \mathbf{g} + \nabla \cdot \tau + \mathbf{F}`

Using the relation :math:`\nabla p = \gamma R_d \pi \nabla \Theta,` where the Exner function :math:`\pi = (p/p_0)^\kappa` with :math:`\kappa = R_d / c_p,`
we can re-write the unapproximated momentum equations

.. math::

  \frac{\partial (\rho u)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} u) - \gamma R_d \pi \frac{\partial \Theta^\prime}{\partial x} + F_u

  \frac{\partial (\rho v)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} v) - \gamma R_d \pi \frac{\partial \Theta^\prime}{\partial y} + F_v

  \frac{\partial (\rho w)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} w) - \gamma R_d \pi \nabla \Theta \frac{\partial \Theta^\prime}{\partial z}
                                                                              - g \overline{\rho} \frac{\pi^\prime}{\overline{\pi}} + F_w
