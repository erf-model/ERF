
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _DryEquations:

Anelastic Equations (Dry)
=============================

ERF can be run in two different modes: in the first, ERF solves the fully compressible fluid equations,
in the second, ERF solves a modified set of equations which approximates the density field with the
hydrostatic density and imposes the anelastic constraint on the velocity field.

In anelastic mode, in the absence of moisture, ERF solves the following partial differential equations
expressing conservation of momentum, potential temperature, and scalars, as well the anelastic constraint
on the velocity.

This option does not currently support terrain-fitted coordinates.

.. math::
  \frac{\partial (\rho_0 \mathbf{u})}{\partial t} &= - \nabla \cdot (\rho_0 \mathbf{u} \mathbf{u}) - \nabla p^\prime
        + \delta_{i,3}\mathbf{B} - \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho_0 \theta)}{\partial t} &= - \nabla \cdot (\rho_0 \mathbf{u} \theta) + \nabla \cdot ( \rho_0 \alpha_{T}\ \nabla \theta) + F_{\rho_0 \theta},

  \frac{\partial (\rho_0 C)}{\partial t} &= - \nabla \cdot (\rho_0 \mathbf{u} C) + \nabla \cdot (\rho_0 \alpha_{C}\ \nabla C)

and

.. math::
  \nabla \cdot \mathbf{u} = 0

where

- :math:`\tau` is the viscous stress tensor,

  .. math::
     \tau_{ij} = -2\mu \sigma_{ij},

with :math:`\sigma_{ij} = S_{ij} -D_{ij}` being the deviatoric part of the strain rate, and

.. math::
   S_{ij} = \frac{1}{2} \left(  \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}   \right), \hspace{24pt}
   D_{ij} = \frac{1}{3}  S_{kk} \delta_{ij} = \frac{1}{3} (\nabla \cdot \mathbf{u}) \delta_{ij},

- :math:`\mathbf{F}` and :math:`F_{\rho \theta}` are the forcing terms described in :ref:`Forcings`,
- :math:`\mathbf{g} = (0,0,-g)` is the gravity vector,

- the potential temperature :math:`\theta` is defined from temperature :math:`T` and hydrostatic pressure :math:`p_0` as

.. math::

  \theta = T \left( \frac{p_0}{p} \right)^{R_d / c_p}.
