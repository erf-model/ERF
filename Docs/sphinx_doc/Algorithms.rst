
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _Algorithms:


Time Advance
============

ERF uses the three-stage low-storage TVD RK3 scheme of Gottlieb and Shu to advance the solution.

Specifically, for

.. math::

  \frac{d \mathbf{U}}{dt} = f(\mathbf{U})

where :math:`\mathbf{U}` is the solution vector, we solve

.. math::

  \mathbf{U}^{n+1/3} = \mathbf{U}^n + \Delta t f(\mathbf{U}^n)

  \mathbf{U}^{n+2/3} = \frac{3}{4} \mathbf{U}^n + \frac{1}{4} ( \mathbf{U}^{n+1/3} + \Delta t f(\mathbf{U}^{n+1/3}) )

  \mathbf{U}^{n+1} = \frac{1}{3} \mathbf{U}^n + \frac{2}{3} ( \mathbf{U}^{n+2/3} + \Delta t f(\mathbf{U}^{n+2/3}) )

In the code, the time-stepping is implemented in :cpp:`Source/RK3/RK3_driver.cpp`, while
:math:`f` is computed in :cpp:`Source/RK3/RK3_stage.cpp`
