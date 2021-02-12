
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

 
.. _Equations:



Compressible Euler Equations
============================

ERF advances the following set of equations:

.. math::
 
  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}),

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + pI) +\rho \mathbf{g},
 
  \frac{\partial (\rho \theta)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \theta),
 
  \frac{\partial (\rho A)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} A),

The relationship between potential temperature and temperature is given by 

.. math::

  \theta = T (\frac{p_0}{p})^{R_d / c_p}

and we use the following equation of state:

.. math::

  p = \rho R_d T;

which can also be written in terms of :math:`\theta` as

.. math::

  p = (\rho R_d \theta / p_0^{R_d / c_p} )^\gamma

Here :math:`\rho, T, \theta`, and :math:`p` are the density, temperature, potential temperature and pressure, respectively; 
these variables are all defined at cell centers.
:math:`A` is an advected quantity, i.e., a tracer, also defined at cell centers.
:math:`\mathbf{u}` and :math:`(\rho \mathbf{u})` are the velocity and momentum, respectively, 
and are defined on faces.  The gravitational vector is denoted by :math:`\mathbf{g}`.

:math:`R_d` and :math:`c_p` are the gas constant and specific heat capacity for dry air respectively, 
and :math:`\gamma = c_p / (c_p - R_d)` .  :math:`p_0` is a reference value for pressure.

