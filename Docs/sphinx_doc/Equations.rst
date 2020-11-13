
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

 
.. _Equations:



Equations
=========

Conservative system
-------------------

ERF advances the following set of fully compressible equations for the conserved state vector: :math:`\mathbf{U} = (\rho, \rho \mathbf{u}, \rho E, \rho Y_k, \rho A_k, \rho B_k):`

.. math::
 
  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}),

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} + pI) +\rho \mathbf{g},

  \frac{\partial (\rho E)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} E + p \mathbf{u}) + \rho \mathbf{u} \cdot \mathbf{g} 
 
  \frac{\partial \rho A}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} A),

and uses the following equation of state:

  p = \rho R_{air} T;

Here :math:`\rho, \mathbf{u}, T`, and :math:`p` are the density, velocity,
temperature and pressure, respectively, and :math:`E
= e + \mathbf{u} \cdot \mathbf{u} / 2` is the total energy with :math:`e` representing the
internal energy.  Here :math:`\mathbf{g}` is the gravitational vector. 
:math:`A_k` is an advected quantity, i.e., a tracer.  

In the code we carry around :math:`T` and :math:`\rho e` in the
state vector even though they are redundant with the state since they may be derived from the other conserved
quantities.  

Some notes:

* There are ``NADV`` advected quantities, which range from :math:`{\tt
  UFA: UFA+nadv-1}`.  The advected quantities have no effect at all on
  the rest of the solution but can be useful as tracer quantities.

Primitive Forms
---------------

ERF uses the primitive form of the fluid equations, defined in terms of
the state :math:`\mathbf{Q} = (\rho, \mathbf{u}, p, \rho e, A)`, to construct the
interface states that are input to the Riemann problem. 
All of the primitive variables are derived from the conservative state
vector. 

The primitive variable equations for density, velocity, and pressure are:

.. math::
  
  \frac{\partial\rho}{\partial t} &=& -\mathbf{u}\cdot\nabla\rho - \rho\nabla\cdot\mathbf{u}

  \frac{\partial\mathbf{u}}{\partial t} &=& -\mathbf{u}\cdot\nabla\mathbf{u} - \frac{1}{\rho}\nabla p + \mathbf{g}

  \frac{\partial p}{\partial t} &=& -\mathbf{u}\cdot\nabla p - \rho c^2\nabla\cdot\mathbf{u}

The advected quantities appear as:

.. math::
  
  \frac{\partial A}{\partial t} &=& -\mathbf{u}\cdot\nabla A 


