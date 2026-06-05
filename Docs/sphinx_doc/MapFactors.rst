
 .. role:: cpp(code)
    :language: c++

.. _sec:MapFactors:

Map Factors (Dry)
=============================

ERF supports the use of isotropic projections only, e.g. Lambert conformal, polar stereographic and Mercator.
Like WRF, ERF implements the projections using map factors, which are defined as the ratio of the distance
in computational space to the corresponding distance on the earth's surface.

With map factors, the dry governing equations have the following form

.. math::
  \frac{\partial \rho}{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u_H} m^{-1})
                                      - m   \frac{\partial (\rho w m^{-1})}{\partial z}

  \frac{\partial (\rho u)}{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u_H} u m^{-1})
                                          - m   \frac{\partial (\rho           w  u m^{-1})}{\partial z}
                                          - m   \frac{\partial p^\prime}{\partial x}
                                          + \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho v)}{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u_H} v m^{-1})
                                          - m   \frac{\partial (\rho           w  v m^{-1})}{\partial z}
                                          - m   \frac{\partial p^\prime}{\partial y}
                                          + \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho w) }{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u_H} w m^{-1})
                                           - m   \frac{\partial (\rho            w w m^{-1})}{\partial z}
                                           -     \frac{\partial p^\prime}{\partial z}
                                           + \rho^\prime \mathbf{g}
                                           + \nabla \cdot \tau + F_z,

  \frac{\partial (\rho \theta)}{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u_H} \theta m^{-1})
                                               - m   \frac{\partial (\rho w \theta m^{-1})}{\partial z}
                                               + \nabla \cdot ( \rho \alpha_{T}\ m \nabla \theta) + F_{\rho \theta},

  \frac{\partial (\rho C)}{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u_H} C m^{-1})
                                          - m   \frac{\partial (\rho w C m^{-1})}{\partial z}
                                          + \nabla \cdot (\rho \alpha_{C}\ m \nabla C)

where
:math:`u_H` is the velocity in the horizontal (lateral) only,
:math:`\nabla_H` is the horizontal (lateral) gradient,
:math:`m` is the map factor at the appropriate spatial location and :math:`m^{-1} = 1 / m`.

The viscous stress tensor
:math:`\tau`
is modified via the strain rates
:math:`S_{ij}`:

.. math::
   S_{11} &= m^2*\left[ \partial_x (um^{-1}) - \frac{h_\xi}{h_\zeta}\partial_z (um^{-1}) \right]

   S_{22} &= m^2*\left[ \partial_y (vm^{-1}) - \frac{h_\eta}{h_\zeta}\partial_z (vm^{-1})) \right]

   S_{33} &= \frac{1}{h_\zeta}\partial_z w

   S_{12} &= \frac{m^2}{2} * \left[ \partial_y (m^{-1}u) + \partial_x (vm^{-1}) - \frac{h_\eta}{h_\zeta} \partial_z (um^{-1})
                                                                                - \frac{h_\xi}{h_\zeta}\partial_z (vm^{-1}) \right]

   S_{13} &= \frac{1}{2} * \left[ \frac{1}{h_\zeta}\partial_z u + m * \left( \partial_x w - \frac{h_\xi}{h_\zeta} \partial_z w \right) \right]

   S_{23} &= \frac{1}{2} * \left[ \frac{1}{h_\zeta}\partial_z v + m * \left( \partial_y w - \frac{h_\eta}{h_\zeta} \partial_z w \right) \right]

When LES models are used, the cell volume is used to compute the eddy viscosities. The cell volume must be adjusted using the map factors:

.. math::
   cellVol = \frac{\Delta x \Delta y \Delta z}{m_x m_y}
