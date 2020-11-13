
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

 
.. _Equations:



Equations
=========

Conservative system
-------------------

ERF advances the following set of fully compressible equations:

.. math::
 
  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}),

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + pI) +\rho \mathbf{g},

  \frac{\partial (\rho E)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} E + p \mathbf{u}) + \rho \mathbf{u} \cdot \mathbf{g} 
 
  \frac{\partial (\rho A)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} A),

and uses the following equation of state:

.. math::

  p = \rho R_{air} T;

Here :math:`\rho, T`, and :math:`p` are the density, temperature and pressure, respectively, 
and :math:`E = e + \mathbf{u} \cdot \mathbf{u} / 2` is the total energy with :math:`e` representing the
internal energy; these variables are all defined at cell centers.
:math:`A` is an advected quantity, i.e., a tracer, also defined at cell centers
\mathbf{u} is the velocity, which is defined on faces.
:math:`\mathbf{g}` is the gravitational vector. 

In the code we carry :math:`T` and :math:`(\rho e)` in the state vector even though they are 
redundant with the state since they may be derived from the other conserved quantities.  

Primitive Forms
---------------

ERF uses the primitive form of the fluid equations to construct the fluxes used to update
the conserved variables. All of the primitive variables are derived from the conservative state
vector. 

The primitive variable equations for density, velocity, and pressure are:

.. math::
  
  \frac{\partial\rho}{\partial t} &=& -\mathbf{u}\cdot\nabla\rho - \rho\nabla\cdot\mathbf{u}

  \frac{\partial\mathbf{u}}{\partial t} &=& -\mathbf{u}\cdot\nabla\mathbf{u} - \frac{1}{\rho}\nabla p + \mathbf{g}

  \frac{\partial p}{\partial t} &=& -\mathbf{u}\cdot\nabla p - \rho c^2\nabla\cdot\mathbf{u}

  \frac{\partial A}{\partial t} &=& -\mathbf{u}\cdot\nabla A 


