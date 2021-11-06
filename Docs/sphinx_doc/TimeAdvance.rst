
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _Algorithms:


Time Advance
============

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
while the boundary conditions for :math:`\mathbf{U}^{(2)` are at time :math:`t^{n+1/2}`.
